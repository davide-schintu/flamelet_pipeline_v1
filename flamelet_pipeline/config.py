from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib
from typing import Any

import numpy as np


CASE_TYPES = {"slf_ideal", "slf_real_supercritical", "slf_vle_subcritical", "fpv_ideal"}
PROPERTY_BACKENDS = {
    "slf_ideal": "ideal",
    "fpv_ideal": "ideal",
    "slf_real_supercritical": "real_supercritical",
    "slf_vle_subcritical": "vle_subcritical",
}


@dataclass(frozen=True)
class PipelineConfig:
    path: Path
    raw: dict[str, Any]

    @property
    def case_type(self) -> str:
        return get(self.raw, "case", "type")

    @property
    def backend(self) -> str:
        return get_optional(self.raw, PROPERTY_BACKENDS[self.case_type], "properties", "backend")

    @property
    def pressure_values(self) -> np.ndarray:
        return as_array(get(self.raw, "streams", "pressure_values"))

    @property
    def output_h5(self) -> Path:
        return resolve_path(self.path, get(self.raw, "output", "h5"))

    @property
    def csv_output_dir(self) -> Path:
        return resolve_path(self.path, get(self.raw, "output", "csv_output_dir"))


def load_config(path: str | Path) -> PipelineConfig:
    config_path = Path(path).expanduser().resolve()
    with config_path.open("rb") as handle:
        raw = tomllib.load(handle)
    validate(raw)
    return PipelineConfig(path=config_path, raw=raw)


def validate(config: dict[str, Any]) -> None:
    case_type = get(config, "case", "type")
    if case_type not in CASE_TYPES:
        raise ValueError(f"[case].type must be one of {sorted(CASE_TYPES)}; got {case_type!r}.")
    pressures = as_array(get(config, "streams", "pressure_values"))
    if pressures.ndim != 1 or pressures.size == 0:
        raise ValueError("[streams].pressure_values must be a non-empty 1D list.")
    if np.any(pressures <= 0.0):
        raise ValueError("[streams].pressure_values must be strictly positive.")
    if case_type.startswith("slf_") or case_type == "fpv_ideal":
        chi_values = as_array(get(config, "flamelet", "chi_values"))
        if chi_values.ndim != 1 or chi_values.size == 0:
            raise ValueError("[flamelet].chi_values must be a non-empty 1D list.")
    if case_type == "fpv_ideal":
        get_optional(config, ["CO2", "H2O"], "fpv", "progress_species")
    get(config, "output", "h5")


def get(config: dict[str, Any], *keys: str) -> Any:
    value: Any = config
    for key in keys:
        value = value[key]
    return value


def get_optional(config: dict[str, Any], default: Any, *keys: str) -> Any:
    value: Any = config
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def as_array(values: Any) -> np.ndarray:
    return np.asarray(values, dtype=float)


def build_grid(section: dict[str, Any], name: str = "grid") -> np.ndarray:
    if "values" in section:
        values = as_array(section["values"])
    else:
        mode = section.get("mode", "logspace_with_zero")
        points = int(section["points"])
        if points < 2:
            raise ValueError(f"{name}.points must be at least 2.")
        if mode == "linspace":
            values = np.linspace(float(section.get("start", 0.0)), float(section.get("stop", 1.0)), points)
        elif mode == "logspace_with_zero":
            start = float(section["start"])
            if start <= 0.0:
                raise ValueError(f"{name}.start must be positive for logspace_with_zero.")
            values = np.concatenate(([0.0], np.logspace(np.log10(start), 0.0, points - 1)))
        else:
            raise ValueError(f"{name}.mode must be 'linspace', 'logspace_with_zero', or explicit values.")
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"{name} must be a non-empty 1D grid.")
    return values


def resolve_path(config_path: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()
