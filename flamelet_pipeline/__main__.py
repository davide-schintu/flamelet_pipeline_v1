from __future__ import annotations

import argparse

from .config import load_config
from .convolution import export_csv
from .runner import run_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description="Unified flamelet pipeline.")
    parser.add_argument("config", help="TOML configuration path.")
    parser.add_argument("--stage", choices=("full", "csv"), default="full")
    args = parser.parse_args()
    if args.stage == "csv":
        export_csv(load_config(args.config))
        return
    run_pipeline(args.config)


if __name__ == "__main__":
    main()
