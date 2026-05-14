from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path
import sys

import numpy as np
from spitfire import ChemicalMechanismSpec, Flamelet, FlameletSpec, Library, Dimension, get_ct_solution_array

from .config import PipelineConfig, build_grid, get, get_optional, resolve_path


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
    higher_requested = requested_chi_sorted[requested_chi_sorted > chi_burn * (1.0 + 1.0e-8)]
    chi_ext = float(higher_requested[0]) if higher_requested.size else chi_burn * float(get_optional(config.raw, 1.25, "fpv", "extinction_chi_factor"))
    chi_ext = float(get_optional(config.raw, chi_ext, "fpv", "extinction_chi"))
    if chi_ext <= chi_burn:
        chi_ext = chi_burn * float(get_optional(config.raw, 1.25, "fpv", "extinction_chi_factor"))
    print(
        f"[fpv] P={pressure:g} Pa: chi_burn={chi_burn:g}, integrating extinguishing branch at chi_ext={chi_ext:g}",
        flush=True,
    )

    burning_slice = Library.squeeze(burning_library[:, int(np.argmax(burning_chi))])
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
