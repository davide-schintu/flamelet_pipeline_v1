# Unified Flamelet Pipeline

Run from this directory:

```bash
python -m flamelet_pipeline prova.toml
```

CSV-only convolution from an existing canonical HDF5:

```bash
python -m flamelet_pipeline prova.toml --stage csv
```

## Vendored Native Modules

The repository contains the source for the two native helpers used by the pipeline.

Build both native modules with:

```bash
PYTHON_BIN=/Users/davide/miniconda3/envs/spenv/bin/python scripts/build_native.sh
```

### PDF convolution helper

```bash
cd flamelet_pipeline/cpp_pdf
make python PYTHON_BIN=/Users/davide/miniconda3/envs/spenv/bin/python
```

The Python extension is loaded from:

```text
flamelet_pipeline/cpp_pdf/build/
```

### RFThermo helper

```bash
cd flamelet_pipeline/rfthermo
python setup.py build_ext --inplace
mkdir -p build
cp RFThermo*.so build/
```

The pipeline first searches:

```text
flamelet_pipeline/rfthermo/build/
flamelet_pipeline/rfthermo/
```

and falls back to any `RFThermo` already installed in the Python environment.

## TOML Notes

If `[flamelet].spitfire_utils` is omitted, the vendored helper
`flamelet_pipeline/spitfire_utils.py` is used.

If `[pdf].build_dir` is omitted, the vendored convolution module under
`flamelet_pipeline/cpp_pdf/build/` is used.

If `[pdf].species` is omitted, all species in the HDF5 database are exported.
