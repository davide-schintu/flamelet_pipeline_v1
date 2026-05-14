from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path
import sys

import numpy as np
import json
import os
from spitfire import ChemicalMechanismSpec, Flamelet, FlameletSpec, Library, Dimension, get_ct_solution_array

from .config import PipelineConfig, build_grid, get, get_optional, resolve_path
from typing import Tuple, Optional


@dataclass
class PressureFlamelet:
    pressure: float
    z: np.ndarray
    coord_name: str
    coord: np.ndarray
    fields: dict[str, np.ndarray]
    species_names: list[str]
    species_y: np.ndarray
    diagnostics: dict[str, np.ndarray]


def run_flamelets_for_pressure(config: PipelineConfig, pressure: float) -> PressureFlamelet:
    if config.case_type == "fpv_ideal":
        return _run_fpv_ideal(config, pressure)
    return _run_slf(config, pressure)


def _load_spitfire_utils(config: PipelineConfig):
    util_path = get_optional(config.raw, None, "flamelet", "spitfire_utils")
    if util_path is not None:
        util_dir = str(resolve_path(config.path, util_path).parent)
        if util_dir not in sys.path:
            sys.path.insert(0, util_dir)
        import spitfire_utils
    else:
        from . import spitfire_utils

    builder = getattr(
        spitfire_utils,
        "build_adiabatic_slfm_library_controlled",
        getattr(spitfire_utils, "build_adiabatic_slfm_library"),
    )
    return builder, spitfire_utils.interpolate_library


def _make_slf_specs(config: PipelineConfig, pressure: float) -> FlameletSpec:
    mech = ChemicalMechanismSpec(get(config.raw, "mechanism", "yaml"), get(config.raw, "mechanism", "phase"))
    oxy = mech.stream(
        "TPX",
        (
            get(config.raw, "streams", "oxidizer_temperature"),
            pressure,
            f"{get(config.raw, 'mechanism', 'oxidizer_species')}:1",
        ),
    )
    fuel = mech.stream(
        "TPX",
        (
            get(config.raw, "streams", "fuel_temperature"),
            pressure,
            f"{get(config.raw, 'mechanism', 'fuel_species')}:1",
        ),
    )
    return FlameletSpec(
        mech_spec=mech,
        initial_condition=get_optional(config.raw, "equilibrium", "flamelet", "initial_condition"),
        oxy_stream=oxy,
        fuel_stream=fuel,
        grid_points=int(get(config.raw, "flamelet", "grid_points")),
    )


def _run_slf(config: PipelineConfig, pressure: float) -> PressureFlamelet:
    builder, interpolate_library = _load_spitfire_utils(config)
    z = build_grid(get(config.raw, "interpolation"), "interpolation")
    requested_chi = np.asarray(get(config.raw, "flamelet", "chi_values"), dtype=float)
    specs = _make_slf_specs(config, pressure)
    cache_root = get_optional(config.raw, None, "flamelet", "chi_output_dir")
    cache_dir = None
    if cache_root:
        cache_dir = str(resolve_path(config.path, cache_root) / f"P_{pressure:g}")

    print(
        f"[flamelet] P={pressure:g} Pa: solving {requested_chi.size} chi values, "
        f"Z grid={z.size}",
        flush=True,
    )
    library = _call_slf_builder(
        builder,
        specs,
        diss_rate_values=requested_chi,
        diss_rate_ref=get_optional(config.raw, "stoichiometric", "flamelet", "diss_rate_ref"),
        verbose=bool(get_optional(config.raw, True, "execution", "verbose")),
        solver_verbose=bool(get_optional(config.raw, False, "execution", "solver_verbose")),
        transient_args=dict(get_optional(config.raw, {}, "flamelet", "transient_args")),
        chi_output_dir=cache_dir,
        solver=get_optional(config.raw, "compute_steady_state", "flamelet", "solver"),
        tolerance=float(get_optional(config.raw, 1.0e-6, "flamelet", "tolerance")),
        use_psitc=bool(get_optional(config.raw, True, "flamelet", "use_psitc")),
        skip_failed=bool(get_optional(config.raw, False, "flamelet", "skip_failed")),
        retry_failed_with_seed=bool(get_optional(config.raw, True, "flamelet", "retry_failed_with_seed")),
        include_extinguished=bool(get_optional(config.raw, False, "flamelet", "include_extinguished")),
    )
    print(f"[flamelet] P={pressure:g} Pa: adding base Cantera properties from the SLF library", flush=True)
    _add_ideal_gas_properties(library, specs.mech_spec, pressure)
    print(f"[flamelet] P={pressure:g} Pa: interpolating library", flush=True)
    interp = interpolate_library(library, z)
    available_coord = np.asarray(library.dissipation_rate_stoich_values, dtype=float)
    order = np.argsort(available_coord)
    available_coord = available_coord[order]
    fields = {name: np.asarray(value)[:, order] for name, value in interp["fields"].items()}
    species_fields = list(interp["species_fields"])
    species_names = [name.replace("mass fraction ", "") for name in species_fields]
    species_y = np.asarray(interp["species_mass_fractions"])[:, order, :]
    coord, fields, species_y, duplicated = _force_requested_chi_axis(
        pressure,
        available_coord,
        np.sort(requested_chi),
        fields,
        species_y,
    )
    diagnostics = {
        "extinguished": np.zeros(coord.shape, dtype=bool),
        "backend_failed": np.zeros((z.size, coord.size), dtype=bool),
        "converged": np.ones((z.size, coord.size), dtype=bool),
        "chi_duplicated": duplicated,
    }
    return PressureFlamelet(pressure, z, "chi", coord, fields, species_names, species_y, diagnostics)


def _call_slf_builder(builder, flamelet_specs: FlameletSpec, **kwargs):
    signature = inspect.signature(builder)
    accepted = signature.parameters
    call_kwargs = {name: value for name, value in kwargs.items() if name in accepted and value is not None}
    if "compute_verbose" in accepted and "compute_verbose" not in call_kwargs:
        call_kwargs["compute_verbose"] = kwargs.get("solver_verbose", False)
    return builder(flamelet_specs, **call_kwargs)


def _force_requested_chi_axis(
    pressure: float,
    available_chi: np.ndarray,
    requested_chi: np.ndarray,
    fields: dict[str, np.ndarray],
    species_y: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, np.ndarray]:
    if available_chi.size == requested_chi.size and np.allclose(available_chi, requested_chi, rtol=1.0e-6, atol=1.0e-12):
        return requested_chi, fields, species_y, np.zeros(requested_chi.shape, dtype=bool)

    source_indices = []
    duplicated = np.zeros(requested_chi.shape, dtype=bool)
    for i, chi in enumerate(requested_chi):
        matches = np.where(np.isclose(available_chi, chi, rtol=1.0e-6, atol=1.0e-12))[0]
        if matches.size:
            source_indices.append(int(matches[0]))
            continue

        lower = np.where(available_chi < chi)[0]
        if lower.size:
            source = int(lower[-1])
        else:
            source = int(np.argmin(np.abs(available_chi - chi)))
        source_indices.append(source)
        duplicated[i] = True
        print(
            f"[flamelet] P={pressure:g} Pa: chi={chi:g} missing, "
            f"duplicating chi={available_chi[source]:g} and storing it as chi={chi:g}",
            flush=True,
        )

    index = np.asarray(source_indices, dtype=int)
    aligned_fields = {name: values[:, index] for name, values in fields.items()}
    aligned_species = species_y[:, index, :]
    return requested_chi, aligned_fields, aligned_species, duplicated


def _add_ideal_gas_properties(library, mechanism: ChemicalMechanismSpec, pressure: float) -> None:
    temperature = np.asarray(library["temperature"])
    if "pressure" not in library:
        library["pressure"] = np.full_like(temperature, pressure, dtype=float)

    ct_states, lib_shape = get_ct_solution_array(mechanism=mechanism, library=library)
    rho = np.asarray(ct_states.density).reshape(lib_shape)
    cp = np.asarray(ct_states.cp_mass).reshape(lib_shape)
    h = np.asarray(ct_states.enthalpy_mass).reshape(lib_shape)
    mu = np.asarray(ct_states.viscosity).reshape(lib_shape)
    lam = np.asarray(ct_states.thermal_conductivity).reshape(lib_shape)
    mean_w = np.asarray(ct_states.mean_molecular_weight).reshape(lib_shape)
    r_mix = 8.31446261815324e3 / mean_w

    library["density"] = rho
    library["enthalpy"] = h
    library["cp"] = cp
    library["viscosity"] = mu
    library["lambda"] = lam
    library["alpha"] = lam / (rho)
    library["Zc"] = np.ones_like(rho)
    library["psi"] = rho / pressure
    library["Rmix"] = r_mix
    library["W"] = mean_w
    library["dpdrho_T"] = np.ones_like(rho)
    if hasattr(ct_states, "sound_speed"):
        sound_speed = np.asarray(ct_states.sound_speed).reshape(lib_shape)
    else:
        cv = np.asarray(ct_states.cv_mass).reshape(lib_shape)
        sound_speed = np.sqrt(np.maximum(1.0e-300, cp / np.maximum(cv, 1.0e-300) * r_mix * temperature))
    library["sound_speed"] = sound_speed


# -----------------------------------------------------------------------------
# Extinction search utilities
# -----------------------------------------------------------------------------

def _is_flamelet_burning(library: Library, T_oxy: float, T_fuel: float, delta_T: float) -> bool:
    """
    Determine whether a flamelet solution represents a burning state.

    Parameters
    ----------
    library:
        Flamelet solution library returned by spitfire.  It is expected to contain
        a "temperature" field.
    T_oxy, T_fuel:
        Temperatures of the oxidizer and fuel inlet streams, used to define a
        baseline cold temperature.
    delta_T:
        Minimum temperature rise above the cold inlet temperature that qualifies
        as burning.  If the maximum flamelet temperature does not exceed
        max(T_oxy, T_fuel) + delta_T the flamelet is considered extinguished.

    Returns
    -------
    bool
        True if the flamelet is burning, False if extinguished or invalid.
    """
    try:
        T = np.asarray(library["temperature"], dtype=float)
    except Exception:
        return False
    if T.size == 0 or not np.isfinite(T).any():
        return False
    Tmax = float(np.nanmax(T))
    if not np.isfinite(Tmax):
        return False
    return Tmax > max(float(T_oxy), float(T_fuel)) + float(delta_T)


def _solve_single_flamelet(
    builder,
    specs: FlameletSpec,
    chi: float,
    config: PipelineConfig,
    pressure: float,
    include_extinguished: bool = True,
    skip_failed: bool = True,
) -> Optional[Library]:
    """
    Solve one steady flamelet at a specified stoichiometric scalar dissipation rate.

    This function is used only during the FPV extinction search. Therefore it can
    use looser solver tolerances than the final manifold construction.
    """

    tolerance = float(
        get_optional(
            config.raw,
            get_optional(config.raw, 1.0e-6, "flamelet", "tolerance"),
            "fpv",
            "extinction_solver_tolerance",
        )
    )

    use_psitc = bool(
        get_optional(
            config.raw,
            get_optional(config.raw, True, "flamelet", "use_psitc"),
            "fpv",
            "extinction_use_psitc",
        )
    )

    solver = get_optional(
        config.raw,
        get_optional(config.raw, "compute_steady_state", "flamelet", "solver"),
        "fpv",
        "extinction_solver",
    )

    transient_args = dict(
        get_optional(
            config.raw,
            get_optional(config.raw, {}, "flamelet", "transient_args"),
            "fpv",
            "extinction_search_transient_args",
        )
    )

    try:
        lib = _call_slf_builder(
            builder,
            specs,
            diss_rate_values=np.asarray([float(chi)], dtype=float),
            diss_rate_ref=get_optional(config.raw, "stoichiometric", "flamelet", "diss_rate_ref"),
            verbose=bool(get_optional(config.raw, True, "execution", "verbose")),
            solver_verbose=bool(get_optional(config.raw, False, "execution", "solver_verbose")),
            transient_args=transient_args,
            solver=solver,
            tolerance=tolerance,
            use_psitc=use_psitc,
            skip_failed=bool(skip_failed),
            retry_failed_with_seed=bool(get_optional(config.raw, True, "flamelet", "retry_failed_with_seed")),
            include_extinguished=bool(include_extinguished),
        )

        if lib.dissipation_rate_stoich_values.size == 0:
            return None

        return lib

    except Exception as exc:
        print(
            f"[fpv] P={pressure:g} Pa: chi={chi:g} failed during extinction search: {exc}",
            flush=True,
        )
        return None



def _load_extinction_database(path: Optional[str]) -> dict:
    """
    Load a JSON database mapping pressures to extinction chi values.

    Parameters
    ----------
    path:
        Path to the database file.  If None or the file does not exist, an empty
        dictionary is returned.

    Returns
    -------
    dict
        Database content as a Python dictionary mapping stringified pressures
        to dictionaries containing ``chi_burn`` and ``chi_ext``.  Additional fields
        may be present.
    """
    if not path:
        return {}
    try:
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as handle:
                return json.load(handle) or {}
    except Exception:
        pass
    return {}


def _save_extinction_database(path: Optional[str], data: dict) -> None:
    """
    Write the extinction chi database back to disk.

    The directory containing the database will be created if it does not exist.

    Parameters
    ----------
    path:
        Path to the JSON file.
    data:
        Dictionary mapping pressures to extinction information.
    """
    if not path:
        return
    try:
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)
    except Exception:
        # If saving fails, silently ignore – the search can still proceed.
        pass


def _find_extinction_chi(
    config: PipelineConfig,
    pressure: float,
    builder,
    specs: FlameletSpec,
    initial_burn_chi: float,
) -> Tuple[float, float]:
    """
    Search for the extinction scalar dissipation rate at a given pressure.

    If an extinction database exists, the cached chi_ext is used only as the
    initial guess. The code still verifies it and rebuilds the burning/extinct
    bracket before returning.
    """

    db_path = get_optional(config.raw, None, "fpv", "extinction_database")
    db_resolved: Optional[str] = None
    if db_path:
        db_resolved = str(resolve_path(config.path, db_path))

    db = _load_extinction_database(db_resolved)
    pressure_key = f"{pressure:g}"

    T_oxy = float(get(config.raw, "streams", "oxidizer_temperature"))
    T_fuel = float(get(config.raw, "streams", "fuel_temperature"))
    delta_T_min = float(get_optional(config.raw, 200.0, "fpv", "burning_delta_T_min"))

    growth = float(get_optional(config.raw, 10.0, "fpv", "extinction_growth_factor"))
    max_growth_steps = int(get_optional(config.raw, 8, "fpv", "extinction_max_growth_steps"))
    max_bisect_steps = int(get_optional(config.raw, 12, "fpv", "extinction_bisect_steps"))
    rel_tol = float(get_optional(config.raw, 0.05, "fpv", "extinction_rel_tol"))

    requested_chi = np.asarray(get(config.raw, "flamelet", "chi_values"), dtype=float)
    last_toml_chi = float(np.max(requested_chi))

    chi_guess = None

    if pressure_key in db:
        record = db[pressure_key]
        try:
            chi_guess = float(record["chi_ext"])
            print(
                f"[fpv] P={pressure:g} Pa: using cached chi_ext={chi_guess:g} as first guess",
                flush=True,
            )
        except Exception:
            chi_guess = None

    if chi_guess is None:
        chi_guess = float(
            get_optional(
                config.raw,
                initial_burn_chi * growth,
                "fpv",
                "extinction_chi_guess",
            )
        )

    if chi_guess <= initial_burn_chi:
        chi_guess = initial_burn_chi * growth
        print(
            f"[fpv] P={pressure:g} Pa: no extinction cache found, starting from chi_guess={chi_guess:g}",
            flush=True,
        )

    def try_chi(val: float) -> Tuple[bool, Optional[Library]]:
        lib = _solve_single_flamelet(
            builder,
            specs,
            val,
            config,
            pressure,
            include_extinguished=True,
            skip_failed=True,
        )

        if lib is None:
            return False, None

        _add_ideal_gas_properties(lib, specs.mech_spec, pressure)
        is_burning = _is_flamelet_burning(lib, T_oxy, T_fuel, delta_T_min)

        print(
            f"[fpv] P={pressure:g} Pa: chi={val:g} -> "
            f"{'burning' if is_burning else 'extinct'}",
            flush=True,
        )

        return is_burning, lib

    # ------------------------------------------------------------------
    # 1. Bracketing around chi_guess
    # ------------------------------------------------------------------

    if chi_guess <= initial_burn_chi:
        chi_guess = initial_burn_chi * growth

    is_burning_guess, _ = try_chi(chi_guess)

    if is_burning_guess:
        chi_low = chi_guess
        chi_high = chi_guess * growth

        for _ in range(max_growth_steps):
            is_burning_high, _ = try_chi(chi_high)

            if not is_burning_high:
                break

            chi_low = chi_high
            chi_high *= growth

        else:
            raise RuntimeError(
                f"Could not bracket extinction at P={pressure:g} Pa. "
                f"chi={chi_high:g} is still burning after {max_growth_steps} growth steps."
            )

    else:
        # chi_guess failed or was classified as extinguished.
        # Do not retest initial_burn_chi with the single-flamelet solver:
        # initial_burn_chi comes from the already-built burning branch and is trusted.
        chi_low = initial_burn_chi
        chi_high = chi_guess

        if chi_high <= chi_low:
            chi_high = chi_low * growth

    if chi_low <= 0.0 or chi_high <= 0.0:
        raise RuntimeError(
            f"Invalid extinction bracket at P={pressure:g} Pa: "
            f"chi_low={chi_low:g}, chi_high={chi_high:g}"
        )

    print(
        f"[fpv] P={pressure:g} Pa: initial extinction bracket "
        f"[{chi_low:g}, {chi_high:g}]",
        flush=True,
    )

    # ------------------------------------------------------------------
    # 2. Logarithmic bisection
    # ------------------------------------------------------------------

    for _ in range(max_bisect_steps):
        if chi_high / chi_low - 1.0 <= rel_tol:
            break

        chi_mid = float(np.sqrt(chi_low * chi_high))
        is_burning_mid, _ = try_chi(chi_mid)

        if is_burning_mid:
            chi_low = chi_mid
        else:
            chi_high = chi_mid

    chi_burn_final = float(chi_low)
    chi_ext_final = float(chi_high)

    print(
        f"[fpv] P={pressure:g} Pa: final extinction bracket "
        f"chi_burn={chi_burn_final:g}, chi_ext={chi_ext_final:g}, "
        f"ratio={chi_ext_final / chi_burn_final:g}",
        flush=True,
    )

    db[pressure_key] = {
        "pressure": float(pressure),
        "chi_burn": chi_burn_final,
        "chi_ext": chi_ext_final,
        "T_oxy": T_oxy,
        "T_fuel": T_fuel,
        "burning_delta_T_min": delta_T_min,
        "extinction_rel_tol": rel_tol,
    }

    _save_extinction_database(db_resolved, db)

    return chi_burn_final, chi_ext_final


def _run_fpv_ideal(config: PipelineConfig, pressure: float) -> PressureFlamelet:
    slf = _run_fpv_source_family(config, pressure)
    product_species = get_optional(config.raw, ["CO2", "H2O", "CO", "H2"], "fpv", "progress_species")
    c_grid = build_grid(get_optional(config.raw, {"mode": "linspace", "points": 32, "start": 0.0, "stop": 1.0}, "fpv", "c_grid"), "fpv.c_grid")
    indices = [slf.species_names.index(name) for name in product_species if name in slf.species_names]
    if not indices:
        raise ValueError(f"None of fpv.progress_species={product_species!r} exists in the flamelet species.")
    yc_source = np.sum(slf.species_y[:, :, indices], axis=2)
    yc_min = np.nanmin(yc_source, axis=0)
    yc_max = np.nanmax(yc_source, axis=0)
    denom = yc_max - yc_min
    if np.any(denom <= 0.0):
        bad = np.where(denom <= 0.0)[0]
        raise RuntimeError(f"FPV progress variable has zero range for source tables {bad.tolist()}.")
    c_source = (yc_source - yc_min[None, :]) / denom[None, :]
    print(
        f"[fpv] P={pressure:g} Pa: remapping {slf.coord.size} source tables "
        f"to square ZxC grid ({slf.z.size} x {c_grid.size})",
        flush=True,
    )
    print(
        f"[fpv] P={pressure:g} Pa: source Yc max range "
        f"{np.nanmin(yc_max):g} .. {np.nanmax(yc_max):g}",
        flush=True,
    )

    fields = {name: _remap_source_tables_to_c(c_source, values, c_grid) for name, values in slf.fields.items()}
    fields["Yc"] = _remap_source_tables_to_c(c_source, yc_source, c_grid)
    species_y = np.stack(
        [_remap_source_tables_to_c(c_source, slf.species_y[:, :, i], c_grid) for i in range(slf.species_y.shape[2])],
        axis=2,
    )
    diagnostics = {
        "backend_failed": np.zeros((slf.z.size, c_grid.size), dtype=bool),
        "converged": np.ones((slf.z.size, c_grid.size), dtype=bool),
    }
    return PressureFlamelet(pressure, slf.z, "C", c_grid, fields, slf.species_names, species_y, diagnostics)


def _run_fpv_source_family(config: PipelineConfig, pressure: float) -> PressureFlamelet:
    builder, interpolate_library = _load_spitfire_utils(config)
    z = build_grid(get(config.raw, "interpolation"), "interpolation")
    requested_chi = np.asarray(get(config.raw, "flamelet", "chi_values"), dtype=float)
    requested_chi_sorted = np.sort(requested_chi)
    specs = _make_slf_specs(config, pressure)

    print(
        f"[fpv] P={pressure:g} Pa: building burning branch for quench detection "
        f"({requested_chi_sorted.size} requested chi values)",
        flush=True,
    )
    burning_library = _call_slf_builder(
        builder,
        specs,
        diss_rate_values=requested_chi_sorted,
        diss_rate_ref=get_optional(config.raw, "stoichiometric", "flamelet", "diss_rate_ref"),
        verbose=bool(get_optional(config.raw, True, "execution", "verbose")),
        solver_verbose=bool(get_optional(config.raw, False, "execution", "solver_verbose")),
        transient_args=dict(get_optional(config.raw, {}, "flamelet", "transient_args")),
        solver=get_optional(config.raw, "compute_steady_state", "flamelet", "solver"),
        tolerance=float(get_optional(config.raw, 1.0e-6, "flamelet", "tolerance")),
        use_psitc=bool(get_optional(config.raw, True, "flamelet", "use_psitc")),
        skip_failed=bool(get_optional(config.raw, True, "fpv", "skip_failed_during_quench_search")),
        retry_failed_with_seed=bool(get_optional(config.raw, True, "flamelet", "retry_failed_with_seed")),
        include_extinguished=False,
    )
    _add_ideal_gas_properties(burning_library, specs.mech_spec, pressure)

    burning_chi = np.asarray(burning_library.dissipation_rate_stoich_values, dtype=float)
    order = np.argsort(burning_chi)
    burning_chi = burning_chi[order]
    if burning_chi.size == 0:
        raise RuntimeError("FPV quench search produced no burning flamelets.")
    chi_burn = float(burning_chi[-1])

    # Determine chi_ext either via explicit configuration or via automated search.
    chi_ext_override = get_optional(config.raw, None, "fpv", "extinction_chi")
    if chi_ext_override is not None:
        chi_ext = float(chi_ext_override)
        if chi_ext <= chi_burn:
            chi_ext = chi_burn * float(get_optional(config.raw, 1.25, "fpv", "extinction_chi_factor"))
        print(
            f"[fpv] P={pressure:g} Pa: chi_burn={chi_burn:g}, using configured chi_ext={chi_ext:g}",
            flush=True,
        )
    else:
        chi_burn, chi_ext = _find_extinction_chi(config, pressure, builder, specs, chi_burn)
        print(
            f"[fpv] P={pressure:g} Pa: chi_burn={chi_burn:g}, integrating extinguishing branch at chi_ext={chi_ext:g}",
            flush=True,
        )

    burning_slice = Library.squeeze(burning_library[:, int(order[-1])])
    burning_slice.extra_attributes["mech_spec"] = specs.mech_spec
    ext_specs = FlameletSpec(library_slice=burning_slice, stoich_dissipation_rate=chi_ext)
    ext_flamelet = Flamelet(ext_specs)
    unsteady = ext_flamelet.integrate_to_steady(
        steady_tolerance=float(get_optional(config.raw, 1.0e-4, "fpv", "extinction_steady_tolerance")),
        **dict(get_optional(config.raw, {}, "fpv", "extinction_transient_args")),
    )
    _add_ideal_gas_properties(unsteady, specs.mech_spec, pressure)

    interp_burning = interpolate_library(burning_library, z)
    fields = {name: np.asarray(values)[:, order] for name, values in interp_burning["fields"].items()}
    species_fields = list(interp_burning["species_fields"])
    species_names = [name.replace("mass fraction ", "") for name in species_fields]
    species_y = np.asarray(interp_burning["species_mass_fractions"])[:, order, :]
    source_coord = burning_chi.copy()

    unsteady_fields, unsteady_species_y = _extract_unsteady_fpv_slices(unsteady, z, species_fields)
    if unsteady_fields:
        for name, values in unsteady_fields.items():
            if name in fields:
                fields[name] = np.concatenate([fields[name], values], axis=1)
        species_y = np.concatenate([species_y, unsteady_species_y], axis=1)
        n_ext = next(iter(unsteady_fields.values())).shape[1]
        ext_coord = chi_ext * (1.0 + np.arange(1, n_ext + 1, dtype=float) / max(n_ext, 1))
        source_coord = np.concatenate([source_coord, ext_coord])
        print(f"[fpv] P={pressure:g} Pa: appended {n_ext} extinguishing transient slices", flush=True)

    diagnostics = {
        "backend_failed": np.zeros((z.size, source_coord.size), dtype=bool),
        "converged": np.ones((z.size, source_coord.size), dtype=bool),
    }
    return PressureFlamelet(pressure, z, "I", source_coord, fields, species_names, species_y, diagnostics)


def _extract_unsteady_fpv_slices(unsteady: Library, z: np.ndarray, species_fields: list[str]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    if len(unsteady.dims) < 2:
        return {}, np.empty((z.size, 0, len(species_fields)))
    dim_names = [dim.name for dim in unsteady.dims]
    mix_axis = next((i for i, name in enumerate(dim_names) if "mixture_fraction" in name), 1)
    z_old = np.asarray(unsteady.mixture_fraction_values, dtype=float)
    fields: dict[str, np.ndarray] = {}
    for prop in unsteady.props:
        arr = np.asarray(unsteady[prop], dtype=float)
        if arr.ndim != 2:
            continue
        if mix_axis == 1:
            arr = arr.T
        if arr.shape[0] != z_old.size:
            continue
        fields[prop] = np.stack([np.interp(z, z_old, arr[:, i]) for i in range(arr.shape[1])], axis=1)
    if species_fields:
        species = np.stack([fields[name] for name in species_fields if name in fields], axis=2)
    else:
        species = np.empty((z.size, 0, 0))
    return fields, species


def _remap_source_tables_to_c(c_source: np.ndarray, values: np.ndarray, c_grid: np.ndarray) -> np.ndarray:
    out = np.empty((values.shape[0], c_grid.size), dtype=float)
    for iz in range(values.shape[0]):
        c_row = np.asarray(c_source[iz, :], dtype=float)
        v_row = np.asarray(values[iz, :], dtype=float)
        good = np.isfinite(c_row) & np.isfinite(v_row)
        if good.sum() == 0:
            out[iz, :] = np.nan
            continue
        if good.sum() == 1:
            out[iz, :] = v_row[good][0]
            continue
        c_good = np.clip(c_row[good], 0.0, 1.0)
        v_good = v_row[good]
        order = np.argsort(c_good)
        c_sorted = c_good[order]
        v_sorted = v_good[order]
        c_unique, inverse = np.unique(c_sorted, return_inverse=True)
        if c_unique.size == 1:
            out[iz, :] = np.mean(v_sorted)
            continue
        v_unique = np.zeros(c_unique.size)
        counts = np.zeros(c_unique.size)
        for i, group in enumerate(inverse):
            v_unique[group] += v_sorted[i]
            counts[group] += 1.0
        v_unique /= counts
        out[iz, :] = np.interp(c_grid, c_unique, v_unique, left=v_unique[0], right=v_unique[-1])
    return out
