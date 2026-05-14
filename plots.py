#!/usr/bin/env python3
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5file")
    parser.add_argument("--pressure-index", "-p", type=int, default=0)
    parser.add_argument("--field", "-f", default="T")
    parser.add_argument("--max-curves", type=int, default=None)
    parser.add_argument("--out", "-o", default=None)
    args = parser.parse_args()

    with h5py.File(args.h5file, "r") as h5:
        print("Root keys:", list(h5.keys()))
        print("Axes:", list(h5["axes"].keys()))
        print("Data fields:", list(h5["data"].keys()))

        Z = h5["axes"]["Z"][:]

        # Per FPV l'asse dovrebbe chiamarsi C.
        if "C" in h5["axes"]:
            C = h5["axes"]["C"][:]
            coord_name = "C"
        elif "chi" in h5["axes"]:
            C = h5["axes"]["chi"][:]
            coord_name = "chi"
        else:
            raise RuntimeError(f"Could not find C or chi axis. Available axes: {list(h5['axes'].keys())}")

        pressures = h5["axes"]["pressure"][:]

        if args.field not in h5["data"]:
            raise RuntimeError(f"Field {args.field!r} not found. Available fields: {list(h5['data'].keys())}")

        data = h5["data"][args.field][:]

    print(f"Loaded field {args.field}: shape = {data.shape}")
    print(f"Pressure axis: {pressures}")
    print(f"Using pressure index {args.pressure_index}: P = {pressures[args.pressure_index]:g} Pa")
    print(f"Coordinate axis: {coord_name}, size = {C.size}")

    # Expected shape: (nP, nZ, nC)
    T = data[args.pressure_index, :, :]

    if args.max_curves is not None and args.max_curves < C.size:
        idx = np.linspace(0, C.size - 1, args.max_curves, dtype=int)
    else:
        idx = np.arange(C.size)

    plt.figure(figsize=(7, 5))

    for j in idx:
        plt.plot(Z, T[:, j], linewidth=1.0, label=f"{coord_name}={C[j]:.3g}")

    plt.xlabel("Z")
    plt.ylabel(args.field)
    plt.title(f"{args.field}(Z) at P = {pressures[args.pressure_index] / 1e5:.2f} bar")
    plt.grid(True, alpha=0.3)

    if len(idx) <= 15:
        plt.legend(fontsize=8)

    plt.tight_layout()

    if args.out:
        plt.savefig(args.out, dpi=300)
        print(f"Saved figure to {args.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()