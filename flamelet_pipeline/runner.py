from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

from .config import PipelineConfig, get_optional, load_config
from .database import assemble_pressure_result, write_hdf5
from .flamelet import run_flamelets_for_pressure
from .properties import compute_properties
from .convolution import export_csv


def run_pipeline(config_path: str) -> str:
    config = load_config(config_path)
    workers = int(get_optional(config.raw, 1, "execution", "workers"))
    pressures = list(config.pressure_values)
    print(f"[pipeline] case={config.case_type}, pressures={len(pressures)}, workers={workers}", flush=True)
    if workers > 1 and len(pressures) > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(_run_pressure_from_path, [(str(config.path), float(p)) for p in pressures]))
    else:
        results = [_run_pressure(config, float(p)) for p in pressures]
    print(f"[pipeline] writing HDF5: {config.output_h5}", flush=True)
    h5_path = write_hdf5(config, results)
    if bool(get_optional(config.raw, True, "execution", "export_csv")):
        print("[pipeline] exporting CSV", flush=True)
        export_csv(config, h5_path)
    print(f"[pipeline] done: {h5_path}", flush=True)
    return h5_path


def _run_pressure_from_path(args: tuple[str, float]) -> dict[str, object]:
    path, pressure = args
    return _run_pressure(load_config(path), pressure)


def _run_pressure(config: PipelineConfig, pressure: float) -> dict[str, object]:
    print(
        f"[pressure] start P={pressure:g} Pa ({pressure / 1.0e5:g} bar), "
        f"backend={config.backend}",
        flush=True,
    )
    flamelet = run_flamelets_for_pressure(config, pressure)
    print(f"[pressure] properties P={pressure:g} Pa", flush=True)
    properties = compute_properties(config, flamelet)
    print(f"[pressure] done P={pressure:g} Pa, {flamelet.coord_name} points={flamelet.coord.size}", flush=True)
    return assemble_pressure_result(flamelet, properties)
