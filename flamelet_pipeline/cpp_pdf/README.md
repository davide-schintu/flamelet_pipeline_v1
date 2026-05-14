# pdf_integrator (C++ + Python bindings)

This folder contains:
- `pdf_integrator`: standalone C++ HDF5 table integrator
- `pdf_integrator_cpp`: Python extension module (pybind11)

## 1) Reproducible environment (recommended)

```bash
conda env create -f environment.yml
conda activate pdf-integrator
```

## 2) Build the C++ executable (portable path)

```bash
cmake -S . -B build-cmake
cmake --build build-cmake -j
```

Optional tuning flags:
- `-DPDF_ENABLE_OPENMP=OFF` to disable OpenMP explicitly
- `-DPDF_ENABLE_NATIVE_OPT=ON` to enable `-march=native` on your own machine only
- `-DPDF_ENABLE_FAST_MATH=ON` to enable `-ffast-math`

## 3) Build Python module

```bash
make python
```

The module is written to `build/pdf_integrator_cpp<ext>.so`.

## 4) Run smoke test

```bash
python tests/test_integrate_1d.py
```

## 5) How to hand off to a colleague

Share source files and metadata, not local build artifacts:
- `src/`, `include/`, `tests/`
- `CMakeLists.txt`, `Makefile`, `environment.yml`, `README.md`

Do **not** share:
- `build/`, `build-cmake/`
- prebuilt `.so` binaries from your machine
