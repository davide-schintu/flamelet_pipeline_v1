from __future__ import annotations

import numpy as np

from .config import PipelineConfig, get, get_optional, resolve_path
from .flamelet import PressureFlamelet
from .thermo import load_rftermo_mixture


UNIVERSAL_GAS_CONSTANT = 8.31446261815324e3
CANONICAL_FIELDS = ("T", "rho", "psi", "Zc", "h", "cp", "mu", "lam", "alpha", "dpdrho_T", "sound_speed", "Rmix", "W")
SPITFIRE_FIELD_ALIASES = {
    "T": ("temperature",),
    "rho": ("density", "density-reynolds"),
    "h": ("enthalpy",),
    "cp": ("cp", "heat capacity cp"),
    "mu": ("viscosity", "mu"),
    "lam": ("lambda", "thermal conductivity", "lam"),
    "alpha": ("alpha",),
    "psi": ("psi",),
    "Zc": ("Zc",),
    "Rmix": ("Rmix",),
    "W": ("W",),
    "dpdrho_T": ("dpdrho_T",),
    "sound_speed": ("sound_speed",),
}


def compute_properties(config: PipelineConfig, flamelet: PressureFlamelet) -> dict[str, object]:
    if config.backend == "ideal":
        return _ideal_properties(config, flamelet)
    if config.backend == "real_supercritical":
        return _rftermo_properties(config, flamelet, vle=False)
    if config.backend == "vle_subcritical":
        return _rftermo_properties(config, flamelet, vle=True)
    raise ValueError(f"Unsupported properties backend {config.backend!r}.")


def _ideal_properties(config: PipelineConfig, flamelet: PressureFlamelet) -> dict[str, object]:
    data: dict[str, np.ndarray] = {}
    for canonical, aliases in SPITFIRE_FIELD_ALIASES.items():
        source = next((name for name in aliases if name in flamelet.fields), None)
        if source is None:
            raise KeyError(
                f"Ideal backend requires one of {aliases!r} for canonical {canonical!r}; "
                f"available fields are {sorted(flamelet.fields)}."
            )
        data[canonical] = np.asarray(flamelet.fields[source], dtype=float)
    for extra in ("Yc", "PV", "omega_c", "HRR"):
        if extra in flamelet.fields:
            data[extra] = np.asarray(flamelet.fields[extra], dtype=float)
    if flamelet.coord_name == "C":
        data["C"] = np.broadcast_to(flamelet.coord[None, :], data["T"].shape).copy()
    return {"data": data, "species_names": flamelet.species_names, "Y": flamelet.species_y, "diagnostics": flamelet.diagnostics}


def _rftermo_properties(config: PipelineConfig, flamelet: PressureFlamelet, vle: bool) -> dict[str, object]:
    Mixture = load_rftermo_mixture()

    mixture_database = str(resolve_path(config.path, get(config.raw, "mechanism", "mixture_database")))
    print(
        f"[properties] P={flamelet.pressure:g} Pa: backend="
        f"{'vle_subcritical' if vle else 'real_supercritical'}, database={mixture_database}",
        flush=True,
    )
    mix = Mixture(mixture_database)
    reorder = _species_reorder_index(flamelet.species_names, mix.Species_names)
    n_coord = flamelet.coord.size
    n_z = flamelet.z.size if flamelet.coord_name == "chi" else flamelet.coord.size
    shape = (n_z, n_coord) if flamelet.coord_name == "chi" else (flamelet.coord.size, 1)
    data = {name: np.full(shape, np.nan) for name in CANONICAL_FIELDS}
    data["T_flamelet"] = np.asarray(flamelet.fields["temperature"], dtype=float).copy()
    data["T_sol"] = np.full(shape, np.nan)
    if vle:
        data["beta"] = np.full(shape, np.nan)
    y_out = np.zeros(shape + (len(mix.Species_names),))
    backend_failed = np.zeros(shape, dtype=bool)
    vle_converged = np.ones(shape, dtype=bool)
    temperature_mode = get_optional(config.raw, "flamelet", "real_gas", "temperature_mode")
    if temperature_mode not in ("flamelet", "enthalpy"):
        raise ValueError("[real_gas].temperature_mode must be 'flamelet' or 'enthalpy'.")
    h_linear = None
    if temperature_mode == "enthalpy":
        h_linear = _build_real_gas_boundary_enthalpy_profile(config, mix, flamelet.pressure, flamelet.z, vle)

    error_count = 0
    max_logged_errors = int(get_optional(config.raw, 8, "execution", "max_backend_error_logs"))
    for j in range(shape[1]):
        for i in range(shape[0]):
            try:
                y = flamelet.species_y[i, j, reorder] if flamelet.coord_name == "chi" else flamelet.species_y[i, 0, reorder]
                t_guess = float(flamelet.fields["temperature"][i, j])
                mix.setTP(t_guess, flamelet.pressure)
                mix.setY(y.tolist())
                if temperature_mode == "enthalpy":
                    h_target = float(h_linear[i])
                    if vle:
                        res_t = mix.solveTemperatureFromH_VLE(h_target, t_guess, 1.0e-4, 500)
                        if isinstance(res_t, tuple):
                            t_sol, is_converged = float(res_t[0]), bool(res_t[1])
                        else:
                            t_sol, is_converged = float(res_t), True
                        if not is_converged or not np.isfinite(t_sol):
                            vle_converged[i, j] = False
                            continue
                    else:
                        t_sol = mix.solveTemperatureFromH(h_target, t_guess, 1.0e-4, 500)
                    mix.setT(t_sol)
                else:
                    t_sol = t_guess
                values = mix.calculateProperties_VLE() if vle else mix.calculateProperties()
                zc, rho, h, cp, mu, lam, alpha, dpdrho_t, sound_speed = values[:9]
                beta = values[9] if vle and len(values) > 9 else np.nan
                vals = np.asarray(values[:10] if vle else values[:9], dtype=float)
                if not np.all(np.isfinite(vals)):
                    if vle:
                        vle_converged[i, j] = False
                    continue
                r_mix = flamelet.pressure / (zc * rho * t_sol)
                data["T"][i, j] = t_sol
                data["T_sol"][i, j] = t_sol
                data["Zc"][i, j] = zc
                data["rho"][i, j] = rho
                data["psi"][i, j] = rho / flamelet.pressure
                data["h"][i, j] = h
                data["cp"][i, j] = cp
                data["mu"][i, j] = mu
                data["lam"][i, j] = lam
                data["alpha"][i, j] = alpha
                data["dpdrho_T"][i, j] = dpdrho_t
                data["sound_speed"][i, j] = sound_speed
                data["Rmix"][i, j] = r_mix
                data["W"][i, j] = UNIVERSAL_GAS_CONSTANT / r_mix
                if vle:
                    data["beta"][i, j] = beta
                    vle_converged[i, j] = True
                y_out[i, j, :] = y
            except Exception as exc:
                backend_failed[i, j] = True
                vle_converged[i, j] = False
                error_count += 1
                if error_count <= max_logged_errors:
                    print(
                        f"[properties] P={flamelet.pressure:g} Pa: backend failed at "
                        f"iZ={i}, i{flamelet.coord_name}={j}, "
                        f"Z={flamelet.z[i] if i < flamelet.z.size else i:g}, "
                        f"{flamelet.coord_name}={flamelet.coord[j] if j < flamelet.coord.size else j:g}, "
                        f"T_guess={locals().get('t_guess', np.nan)!r}: "
                        f"{type(exc).__name__}: {exc}",
                        flush=True,
                    )

    if not np.isfinite(data["T"]).any():
        raise RuntimeError(
            f"RFThermo backend produced no finite points for P={flamelet.pressure:g} Pa. "
            "Check [mechanism].mixture_database, species ordering, and the logged backend errors above."
        )

    if vle:
        _fill_internal_nans(data, exclude={"T_flamelet"}, valid=vle_converged)
        data["beta"] = np.clip(data["beta"], 0.0, 1.0)
    diagnostics = dict(flamelet.diagnostics)
    diagnostics["backend_failed"] = backend_failed
    if vle:
        diagnostics["vle_converged"] = vle_converged
    return {"data": data, "species_names": list(mix.Species_names), "Y": y_out, "diagnostics": diagnostics}


def _build_real_gas_boundary_enthalpy_profile(config: PipelineConfig, mix, pressure: float, z: np.ndarray, vle: bool) -> np.ndarray:
    species_names = list(mix.Species_names)
    fuel_species = get(config.raw, "mechanism", "fuel_species")
    oxidizer_species = get(config.raw, "mechanism", "oxidizer_species")
    fuel_temperature = float(get(config.raw, "streams", "fuel_temperature"))
    oxidizer_temperature = float(get(config.raw, "streams", "oxidizer_temperature"))
    h_fuel = _pure_stream_enthalpy(mix, species_names, fuel_species, fuel_temperature, pressure, vle)
    h_oxidizer = _pure_stream_enthalpy(mix, species_names, oxidizer_species, oxidizer_temperature, pressure, vle)
    print(
        f"[properties] enthalpy mode: h({oxidizer_species}@{oxidizer_temperature:g} K)={h_oxidizer:g}, "
        f"h({fuel_species}@{fuel_temperature:g} K)={h_fuel:g}",
        flush=True,
    )
    return z * h_fuel + (1.0 - z) * h_oxidizer


def _pure_stream_enthalpy(mix, species_names: list[str], species_name: str, temperature: float, pressure: float, vle: bool) -> float:
    index = next((i for i, name in enumerate(species_names) if name.upper() == species_name.upper()), None)
    if index is None:
        raise ValueError(f"Species {species_name!r} not found in RFThermo database species: {species_names}")
    y = np.zeros(len(species_names))
    y[index] = 1.0
    mix.setTP(temperature, pressure)
    mix.setY(y.tolist())
    values = mix.calculateProperties_VLE() if vle else mix.calculateProperties()
    return float(values[2])


def _species_reorder_index(source_names: list[str], target_names: list[str]) -> list[int]:
    upper = {name.upper(): i for i, name in enumerate(source_names)}
    return [upper[name.upper()] for name in target_names]


def _fill_internal_nans(data: dict[str, np.ndarray], exclude: set[str], valid: np.ndarray | None = None) -> None:
    try:
        from scipy.interpolate import PchipInterpolator
    except Exception:
        PchipInterpolator = None
    for name, values in data.items():
        if name in exclude or values.ndim != 2:
            continue
        x = np.arange(values.shape[0], dtype=float)
        for j in range(values.shape[1]):
            col = values[:, j]
            good = np.isfinite(col)
            if valid is not None:
                good &= np.asarray(valid[:, j], dtype=bool)
            if good.sum() < 2:
                continue
            first = int(np.argmax(good))
            last = len(good) - int(np.argmax(good[::-1]))
            gap = ~good[first:last]
            if not np.any(gap):
                continue
            if PchipInterpolator is not None and good.sum() >= 3:
                col[first:last][gap] = PchipInterpolator(x[good], col[good])(x[first:last][gap])
            else:
                col[first:last][gap] = np.interp(x[first:last][gap], x[good], col[good])
