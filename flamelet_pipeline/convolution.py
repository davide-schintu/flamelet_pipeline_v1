from __future__ import annotations

import os
from pathlib import Path
import sys

import h5py
import numpy as np

from .config import PipelineConfig, as_array, get, get_optional, resolve_path


DIRECT_PDF_FIELDS = {"Zc", "mu", "lam", "alpha", "sound_speed", "dpdrho_T", "Rmix"}
HARMONIC_FIELDS = {"rho", "psi"}
FAVRE_FIELDS = {"T", "h", "cp"}
PDF_AVERAGING_MODES = {"opensmoke", "reynolds"}


def export_csv(config: PipelineConfig, db_path: str | os.PathLike[str] | None = None) -> None:
    pdf = _load_pdf_integrator(config)
    averaging = get_optional(config.raw, "opensmoke", "pdf", "averaging")
    if averaging not in PDF_AVERAGING_MODES:
        raise ValueError("[pdf].averaging must be 'opensmoke' or 'reynolds'.")
    _require_pdf_functions(pdf)
    h5_path = Path(db_path) if db_path is not None else config.output_h5
    s_vector = as_array(get(config.raw, "pdf", "s_vector"))
    configured_species = get_optional(config.raw, None, "pdf", "species")
    output_dir = config.csv_output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "r") as handle:
        pressure = handle["axes/pressure"][:]
        z = handle["axes/Z"][:]
        coord_name = handle.attrs.get("coord_name", "chi")
        coord = handle[f"axes/{coord_name}"][:]
        species_names = [name.decode() if isinstance(name, bytes) else str(name) for name in handle["species/names"][:]]
        requested_species = list(species_names if configured_species is None else configured_species)
        species_index = {name.upper(): i for i, name in enumerate(species_names)}
        csv_fields = list(get(config.raw, "output", "csv_fields")) + requested_species
        header = ",".join(csv_fields)
        for ip, pressure_value in enumerate(pressure):
            p_dir = output_dir / f"P_{pressure_value / 1.0e5:g}bar"
            print(
                f"[csv] P={pressure_value:g} Pa ({pressure_value / 1.0e5:g} bar): "
                f"{coord_name} points={coord.size}, ZVar cases={s_vector.size}",
                flush=True,
            )
            for ic, coord_value in enumerate(coord):
                if coord_name == "chi":
                    coord_dir = f"Chi_{coord_value:g}"
                elif coord_name == "C":
                    coord_dir = f"C_{coord_value:g}"
                else:
                    coord_dir = f"{coord_name}_{coord_value:g}"
                c_dir = p_dir / coord_dir
                print(f"[csv] P={pressure_value:g} Pa: {coord_name}={coord_value:g}", flush=True)
                c_dir.mkdir(parents=True, exist_ok=True)
                rho_base = handle["data/rho"][ip, :, ic]
                rho3 = rho_base[None, None, :]
                field_base = _collect_fields(handle, ip, ic, csv_fields, requested_species, species_index)
                for s_value in s_vector:
                    columns = []
                    rho_int = pdf.integrate_density_opensmoke(z, np.asarray([s_value]), rho3)[0, 0, :, 0]
                    harmonic_cache = {"rho": rho_int}
                    for field in csv_fields:
                        if field == "Z":
                            data = z
                        elif field == coord_name:
                            data = np.full_like(z, coord_value, dtype=float)
                        elif field in HARMONIC_FIELDS:
                            if field not in harmonic_cache:
                                harmonic_cache[field] = pdf.integrate_density_opensmoke(
                                    z,
                                    np.asarray([s_value]),
                                    field_base[field][None, None, :],
                                )[0, 0, :, 0]
                            data = harmonic_cache[field]
                        elif averaging == "opensmoke" and (field in FAVRE_FIELDS or field in requested_species):
                            data = pdf.integrate_favre(z, np.asarray([s_value]), rho3, field_base[field][None, None, :])[0, 0, :, 0]
                        elif averaging == "reynolds" or field in DIRECT_PDF_FIELDS:
                            data = _integrate_reynolds(pdf, z, field_base[field], s_value)
                        else:
                            raise ValueError(f"No convolution rule for field {field!r}.")
                        columns.append(data)
                    np.savetxt(c_dir / f"ZVar_{s_value:g}.csv", np.column_stack(columns), delimiter=",", header=header, comments="")


def _load_pdf_integrator(config: PipelineConfig):
    build_dir_cfg = get_optional(config.raw, None, "pdf", "build_dir")
    if build_dir_cfg is None:
        build_dir = Path(__file__).resolve().parent / "cpp_pdf" / "build"
    else:
        build_dir = resolve_path(config.path, build_dir_cfg)
    if str(build_dir) not in sys.path:
        sys.path.insert(0, str(build_dir))
    import pdf_integrator_cpp

    return pdf_integrator_cpp


def _require_pdf_functions(pdf) -> None:
    missing = [name for name in ("integrate_density_opensmoke", "integrate_favre") if not hasattr(pdf, name)]
    if missing:
        raise RuntimeError(f"pdf_integrator_cpp is missing required OpenSMOKE functions: {', '.join(missing)}")


def _integrate_reynolds(pdf, z: np.ndarray, field: np.ndarray, s_value: float) -> np.ndarray:
    zvar = s_value * (z * (1.0 - z))
    return pdf.integrate_1d(z, field, zvar)


def _collect_fields(handle, ip: int, ic: int, fields: list[str], requested_species: list[str], species_index: dict[str, int]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for field in fields:
        if field == "Z":
            continue
        if field in requested_species:
            idx = species_index.get(field.upper())
            out[field] = handle["species/Y"][ip, :, ic, idx] if idx is not None else np.zeros(handle["axes/Z"].shape)
        elif f"data/{field}" in handle:
            out[field] = handle[f"data/{field}"][ip, :, ic]
        else:
            raise ValueError(f"Field {field!r} not found in canonical HDF5.")
    return out
