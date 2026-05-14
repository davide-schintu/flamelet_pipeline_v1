from __future__ import annotations

import h5py
import numpy as np

from .config import PipelineConfig


def write_hdf5(config: PipelineConfig, pressure_results: list[dict[str, object]]) -> str:
    pressure_results = sorted(pressure_results, key=lambda item: float(item["pressure"]))
    pressures = np.asarray([item["pressure"] for item in pressure_results], dtype=float)
    first = pressure_results[0]
    z = np.asarray(first["z"], dtype=float)
    coord_name = str(first["coord_name"])
    coord = np.asarray(first["coord"], dtype=float)
    out_path = config.output_h5
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_path, "w") as handle:
        axes = handle.create_group("axes")
        axes.create_dataset("pressure", data=pressures)
        axes.create_dataset("Z", data=z)
        axes.create_dataset(coord_name, data=coord)

        data_group = handle.create_group("data")
        for field in first["data"]:
            stack = np.stack([np.asarray(item["data"][field]) for item in pressure_results], axis=0)
            data_group.create_dataset(field, data=stack)

        species_group = handle.create_group("species")
        species_names = list(first["species_names"])
        species_group.create_dataset("names", data=np.asarray(species_names, dtype="S"))
        species_group.create_dataset("Y", data=np.stack([np.asarray(item["Y"]) for item in pressure_results], axis=0))

        diag_group = handle.create_group("diagnostics")
        diagnostic_names = set().union(*(item["diagnostics"].keys() for item in pressure_results))
        for name in sorted(diagnostic_names):
            stack = np.stack([np.asarray(item["diagnostics"].get(name, np.zeros_like(first["data"]["T"], dtype=bool))) for item in pressure_results], axis=0)
            diag_group.create_dataset(name, data=stack)

        handle.attrs["case_type"] = config.case_type
        handle.attrs["coord_name"] = coord_name
    return str(out_path)


def assemble_pressure_result(flamelet, properties: dict[str, object]) -> dict[str, object]:
    return {
        "pressure": flamelet.pressure,
        "z": flamelet.z,
        "coord_name": flamelet.coord_name,
        "coord": flamelet.coord,
        "data": properties["data"],
        "species_names": properties["species_names"],
        "Y": properties["Y"],
        "diagnostics": properties["diagnostics"],
    }
