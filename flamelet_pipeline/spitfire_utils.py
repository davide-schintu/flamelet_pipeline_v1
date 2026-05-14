from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.interpolate import interp1d

from spitfire.chemistry.flamelet import Flamelet, FlameletSpec
from spitfire.chemistry.library import Dimension, Library


MIXTURE_FRACTION_NAME = "mixture_fraction"
DISSIPATION_RATE_NAME = "dissipation_rate"
STOICH_SUFFIX = "_stoich"


def build_adiabatic_slfm_library_controlled(
    flamelet_specs,
    diss_rate_values=np.logspace(-3, 2, 16),
    diss_rate_ref="stoichiometric",
    verbose=True,
    solver_verbose=False,
    include_extinguished=False,
    diss_rate_log_scaled=True,
    tolerance=1.0e-6,
    solver="compute_steady_state",
    use_psitc=True,
    newton_args=None,
    psitc_args=None,
    transient_args=None,
    skip_failed=False,
    retry_failed_with_seed=False,
    return_intermediates=False,
):
    """Build an adiabatic SLFM library with explicit solver controls.

    This mirrors Spitfire's build_adiabatic_slfm_library loop, but exposes the
    solver arguments hidden by the stock wrapper.
    """
    if isinstance(flamelet_specs, dict):
        flamelet_specs = FlameletSpec(**flamelet_specs)

    solver = solver.lower()
    valid_solvers = {"compute_steady_state", "newton", "psitc", "transient"}
    if solver not in valid_solvers:
        raise ValueError(f"solver must be one of {sorted(valid_solvers)}, got {solver!r}")

    if diss_rate_ref not in {"stoichiometric", "maximum"}:
        raise ValueError("diss_rate_ref must be 'stoichiometric' or 'maximum'")

    flamelet_specs.initial_condition = "equilibrium"
    if diss_rate_ref == "maximum":
        flamelet_specs.max_dissipation_rate = 0.0
        suffix = "_max"
    else:
        flamelet_specs.stoich_dissipation_rate = 0.0
        suffix = STOICH_SUFFIX

    first_flamelet = Flamelet(flamelet_specs)
    table_dict = {}
    chi_st_values = []
    pending_retry_values = []
    diss_rate_values = np.atleast_1d(np.asarray(diss_rate_values, dtype=float))

    if verbose:
        mech = flamelet_specs.mech_spec
        print("-" * 82)
        print("building controlled adiabatic SLFM library")
        print("-" * 82)
        print(f"- mechanism: {mech.mech_file_path}")
        print(f"- {mech.n_species} species, {mech.n_reactions} reactions")
        print(
            "- stoichiometric mixture fraction: "
            f"{mech.stoich_mixture_fraction(flamelet_specs.fuel_stream, flamelet_specs.oxy_stream):.3f}"
        )
        print("-" * 82)

    cput00 = perf_counter()
    for idx, chi_input in enumerate(diss_rate_values):
        if _has_chi_input(table_dict, chi_input):
            continue

        result = _solve_one_chi(
            idx=idx,
            count=diss_rate_values.size,
            chi_input=chi_input,
            suffix=suffix,
            flamelet_specs=flamelet_specs,
            diss_rate_ref=diss_rate_ref,
            solver=solver,
            tolerance=tolerance,
            solver_verbose=solver_verbose,
            use_psitc=use_psitc,
            newton_args=newton_args,
            psitc_args=psitc_args,
            transient_args=transient_args,
            include_extinguished=include_extinguished,
            verbose=verbose,
        )

        if not result["success"]:
            if retry_failed_with_seed:
                pending_retry_values.append(chi_input)
                continue
            if not skip_failed:
                if result["exception"] is not None:
                    raise result["exception"]
                break
            continue

        _add_solution(
            result,
            table_dict,
            chi_st_values,
            flamelet_specs,
            return_intermediates,
        )

        if retry_failed_with_seed and pending_retry_values:
            retry_values = pending_retry_values[::-1]
            pending_retry_values = []
            for retry_chi in retry_values:
                retry_result = _solve_one_chi(
                    idx=None,
                    count=diss_rate_values.size,
                    chi_input=retry_chi,
                    suffix=suffix,
                    flamelet_specs=flamelet_specs,
                    diss_rate_ref=diss_rate_ref,
                    solver=solver,
                    tolerance=tolerance,
                    solver_verbose=solver_verbose,
                    use_psitc=use_psitc,
                    newton_args=newton_args,
                    psitc_args=psitc_args,
                    transient_args=transient_args,
                    include_extinguished=include_extinguished,
                    verbose=verbose,
                    label="retry",
                )
                if retry_result["success"]:
                    _add_solution(
                        retry_result,
                        table_dict,
                        chi_st_values,
                        flamelet_specs,
                        return_intermediates,
                    )
                elif not skip_failed:
                    if retry_result["exception"] is not None:
                        raise retry_result["exception"]
                    break

    if verbose:
        print("-" * 82)
        print(f"library built in {perf_counter() - cput00:6.2f} s")
        print("-" * 82, flush=True)

    if return_intermediates:
        return table_dict, first_flamelet.mixfrac_grid, np.asarray(chi_st_values)

    if not chi_st_values:
        raise RuntimeError(
            "No flamelets were added to the library. Every flamelet either failed or was classified as extinguished."
        )

    chi_st_values = np.asarray(sorted(chi_st_values), dtype=float)
    z_dim = Dimension(MIXTURE_FRACTION_NAME, first_flamelet.mixfrac_grid)
    chi_dim = Dimension(DISSIPATION_RATE_NAME + STOICH_SUFFIX, chi_st_values, diss_rate_log_scaled)
    output_library = Library(z_dim, chi_dim)
    output_library.extra_attributes["mech_spec"] = flamelet_specs.mech_spec

    for quantity in table_dict[chi_st_values[-1]]:
        if quantity.startswith("_"):
            continue
        output_library[quantity] = output_library.get_empty_dataset()
        for ix, chi_st in enumerate(chi_st_values):
            output_library[quantity][:, ix] = table_dict[chi_st][quantity]

    return output_library


def _solve_one_chi(
    idx,
    count,
    chi_input,
    suffix,
    flamelet_specs,
    diss_rate_ref,
    solver,
    tolerance,
    solver_verbose,
    use_psitc,
    newton_args,
    psitc_args,
    transient_args,
    include_extinguished,
    verbose,
    label=None,
):
    if diss_rate_ref == "maximum":
        flamelet_specs.max_dissipation_rate = chi_input
    else:
        flamelet_specs.stoich_dissipation_rate = chi_input

    flamelet = Flamelet(flamelet_specs)
    if verbose:
        if idx is None:
            prefix = f"{label:>8}" if label else "   retry"
        else:
            prefix = f"{idx + 1:4}/{count:4}"
        print(
            f"{prefix} (chi{suffix} = {chi_input:8.1e} 1/s) ",
            end="",
            flush=True,
        )

    cput0 = perf_counter()
    try:
        x_library = _solve_flamelet(
            flamelet,
            solver=solver,
            tolerance=tolerance,
            solver_verbose=solver_verbose,
            use_psitc=use_psitc,
            newton_args=newton_args,
            psitc_args=psitc_args,
            transient_args=transient_args,
        )
    except Exception as exc:
        if verbose:
            print(f" failed ({type(exc).__name__}: {exc}), waiting for a better seed.")
        return {"success": False, "exception": exc}

    dcput = perf_counter() - cput0
    temperature_rise = np.max(flamelet.current_temperature) - np.max(flamelet.linear_temperature)
    if temperature_rise < 10.0 and not include_extinguished:
        if verbose:
            print(
                " extinction detected "
                f"(max temperature rise = {temperature_rise:.2f} K), waiting for a better seed."
            )
        return {"success": False, "exception": None}

    if verbose:
        iteration_count = getattr(flamelet, "_iteration_count", None)
        iter_msg = f", iterations = {iteration_count}" if iteration_count is not None else ""
        print(f" converged in {dcput:6.2f} s{iter_msg}, T_max = {np.max(flamelet.current_temperature):6.1f}")

    return {"success": True, "flamelet": flamelet, "library": x_library, "chi_input": chi_input}


def _add_solution(result, table_dict, chi_st_values, flamelet_specs, return_intermediates):
    flamelet = result["flamelet"]
    x_library = result["library"]
    z_st = flamelet.mechanism.stoich_mixture_fraction(flamelet.fuel_stream, flamelet.oxy_stream)
    chi_st = flamelet._compute_dissipation_rate(
        np.array([z_st]),
        flamelet._max_dissipation_rate,
        flamelet._dissipation_rate_form,
    )[0]

    if chi_st not in table_dict:
        chi_st_values.append(chi_st)
    table_dict[chi_st] = {}
    table_dict[chi_st]["_chi_input"] = result["chi_input"]
    for prop in x_library.props:
        table_dict[chi_st][prop] = x_library[prop].ravel()

    flamelet_specs.initial_condition = flamelet.current_interior_state
    if return_intermediates:
        table_dict[chi_st]["adiabatic_state"] = np.copy(flamelet.current_interior_state)


def _has_chi_input(table_dict, chi_input):
    for data in table_dict.values():
        if "_chi_input" in data and np.isclose(data["_chi_input"], chi_input):
            return True
    return False


def _coerce_flamelet_specs(flamelet_specs):
    if isinstance(flamelet_specs, FlameletSpec):
        return flamelet_specs
    if isinstance(flamelet_specs, dict):
        return FlameletSpec(**flamelet_specs)
    raise TypeError("flamelet_specs must be a FlameletSpec or a dict of FlameletSpec kwargs.")


def _coerce_dissipation_rate_values(diss_rate_values):
    values = np.atleast_1d(np.asarray(diss_rate_values, dtype=float))
    if values.ndim != 1 or values.size == 0:
        raise ValueError("diss_rate_values must be a non-empty 1D array-like.")
    return values


def _chi_filename(index, requested_chi):
    safe_value = f"{requested_chi:.6e}".replace("+", "")
    return f"chi_{index:04d}_{safe_value}.npz"


def save_library_slice(output_dir, index, requested_chi, actual_chi, x_library, mixfrac_grid):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "requested_chi": np.array([requested_chi], dtype=float),
        "chi_st": np.array([actual_chi], dtype=float),
        "mixture_fraction": np.asarray(mixfrac_grid, dtype=float),
    }
    for prop in x_library.props:
        payload[prop] = np.asarray(x_library[prop]).ravel()

    np.savez_compressed(output_dir / _chi_filename(index, requested_chi), **payload)


def load_library_slice(output_dir, index, requested_chi):
    slice_path = Path(output_dir) / _chi_filename(index, requested_chi)
    if not slice_path.exists():
        return None
    with np.load(slice_path) as payload:
        return {key: payload[key].copy() for key in payload.files}


def build_output_library(
    table_dict,
    chi_values,
    mixfrac_grid,
    mech_spec=None,
    flamelet_specs=None,
    diss_rate_log_scaled=True,
    diss_rate_ref="stoichiometric",
):
    chi_values = _coerce_dissipation_rate_values(chi_values)
    suffix = STOICH_SUFFIX if diss_rate_ref == "stoichiometric" else "_max"

    z_dim = Dimension(MIXTURE_FRACTION_NAME, np.asarray(mixfrac_grid, dtype=float))
    chi_dim = Dimension(DISSIPATION_RATE_NAME + suffix, chi_values, diss_rate_log_scaled)
    output_library = Library(z_dim, chi_dim)

    if mech_spec is None and flamelet_specs is not None:
        mech_spec = flamelet_specs.mech_spec
    if mech_spec is not None:
        output_library.extra_attributes["mech_spec"] = mech_spec
    if flamelet_specs is not None:
        output_library.extra_attributes["dissipation_rate_form"] = flamelet_specs.dissipation_rate_form
        output_library.extra_attributes["fuel_stream"] = flamelet_specs.fuel_stream
        output_library.extra_attributes["oxy_stream"] = flamelet_specs.oxy_stream

    first_chi = chi_values[0]
    for quantity in table_dict[first_chi]:
        if quantity.startswith("_"):
            continue
        output_library[quantity] = output_library.get_empty_dataset()
        for ix, chi_value in enumerate(chi_values):
            output_library[quantity][:, ix] = table_dict[chi_value][quantity]

    return output_library


def load_cached_library(output_dir, diss_rate_values, flamelet_specs, mixfrac_grid=None, diss_rate_log_scaled=True,
                        diss_rate_ref="stoichiometric"):
    diss_rate_values = _coerce_dissipation_rate_values(diss_rate_values)
    cached_slices = []
    for index, requested_chi in enumerate(diss_rate_values):
        cached_slice = load_library_slice(output_dir, index, requested_chi)
        if cached_slice is None:
            return None
        cached_slices.append(cached_slice)

    if mixfrac_grid is None:
        mixfrac_grid = np.asarray(cached_slices[0]["mixture_fraction"], dtype=float)

    chi_values = [float(cached_slice["chi_st"][0]) for cached_slice in cached_slices]
    table_dict = {}
    for chi_value, cached_slice in zip(chi_values, cached_slices):
        table_dict[chi_value] = {
            key: np.asarray(value).ravel()
            for key, value in cached_slice.items()
            if key not in {"requested_chi", "chi_st", "mixture_fraction"}
        }

    return build_output_library(
        table_dict=table_dict,
        chi_values=chi_values,
        mixfrac_grid=mixfrac_grid,
        flamelet_specs=flamelet_specs,
        diss_rate_log_scaled=diss_rate_log_scaled,
        diss_rate_ref=diss_rate_ref,
    )


def interpolate_library_field(field, z_old, z_new, kind="linear"):
    interpolator = interp1d(
        z_old,
        field,
        axis=0,
        kind=kind,
        bounds_error=False,
        fill_value="extrapolate",
    )
    return interpolator(z_new)


def interpolate_library(library, z_new, field_names=None, species_prefix="mass fraction ", kind="linear"):
    z_old = np.asarray(library.mixture_fraction_values, dtype=float)
    z_new = np.asarray(z_new, dtype=float)

    if field_names is None:
        field_names = list(library.props)

    interpolated = {
        field_name: interpolate_library_field(library[field_name], z_old, z_new, kind=kind)
        for field_name in field_names
    }

    species_fields = [field_name for field_name in field_names if field_name.startswith(species_prefix)]
    species_stack = None
    if species_fields:
        species_stack = np.stack([interpolated[field_name] for field_name in species_fields], axis=2)

    return {
        "z": z_new,
        "fields": interpolated,
        "species_fields": species_fields,
        "species_mass_fractions": species_stack,
    }


def _solve_flamelet(
    flamelet,
    solver,
    tolerance,
    solver_verbose,
    use_psitc,
    newton_args,
    psitc_args,
    transient_args,
):
    if solver == "compute_steady_state":
        return flamelet.compute_steady_state(
            tolerance=tolerance,
            verbose=solver_verbose,
            use_psitc=use_psitc,
            newton_args=newton_args,
            psitc_args=psitc_args,
            transient_args=transient_args,
        )

    if solver == "newton":
        args = {"tolerance": tolerance, "verbose": solver_verbose}
        if newton_args is not None:
            args.update(newton_args)
        library, _, converged = flamelet.steady_solve_newton(**args)
        if not converged:
            raise RuntimeError("Newton steady solve did not converge.")
        return library

    if solver == "psitc":
        args = {"tolerance": tolerance, "verbose": solver_verbose}
        if psitc_args is not None:
            args.update(psitc_args)
        library, _, converged, _ = flamelet.steady_solve_psitc(**args)
        if not converged:
            raise RuntimeError("PsiTC steady solve did not converge.")
        return library

    args = {
        "steady_tolerance": tolerance,
        "write_log": solver_verbose,
        "save_first_and_last_only": True,
    }
    if transient_args is not None:
        args.update(transient_args)
    transient_library = flamelet.integrate_to_steady(**args)

    steady_library = Library(transient_library.dim(MIXTURE_FRACTION_NAME))
    steady_library.extra_attributes["mech_spec"] = flamelet.mechanism
    for prop in transient_library.props:
        steady_library[prop] = transient_library[prop][-1, :].ravel()
    return steady_library


build_adiabatic_slfm_library = build_adiabatic_slfm_library_controlled
