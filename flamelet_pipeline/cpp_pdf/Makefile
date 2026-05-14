# Makefile for PDF Integrator

CXX ?= c++
BASE_CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra
NATIVE_OPT ?= 0
FAST_MATH ?= 0

ifeq ($(NATIVE_OPT),1)
BASE_CXXFLAGS += -march=native
endif
ifeq ($(FAST_MATH),1)
BASE_CXXFLAGS += -ffast-math
endif

PYTHON_BIN ?= python3
PY_PREFIX = $(shell $(PYTHON_BIN) -c "import pathlib,sys; print(pathlib.Path(sys.executable).resolve().parents[1])")

# Override this if HDF5/OpenMP are in a different prefix
CONDA_PREFIX ?= $(PY_PREFIX)

HDF5_INC = $(CONDA_PREFIX)/include
HDF5_LIB = $(CONDA_PREFIX)/lib
OPENMP_INC = $(CONDA_PREFIX)/include
OPENMP_LIB = $(CONDA_PREFIX)/lib

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
OMP_CFLAGS = -Xclang -fopenmp
OMP_LDFLAGS = -lomp
PY_SHARED_LDFLAGS = -undefined dynamic_lookup
OPENMP_LIB_FILE = libomp.dylib
else
OMP_CFLAGS = -fopenmp
OMP_LDFLAGS = -fopenmp
PY_SHARED_LDFLAGS =
OPENMP_LIB_FILE = libgomp.so
endif

INCLUDES = -Iinclude -I$(HDF5_INC) -I$(OPENMP_INC)
LDFLAGS = -L$(HDF5_LIB) -L$(OPENMP_LIB) -lhdf5_cpp -lhdf5 $(OMP_LDFLAGS) -Wl,-rpath,$(HDF5_LIB)

PY_INCLUDES = $(shell $(PYTHON_BIN) -m pybind11 --includes)
PY_CFLAGS = $(shell $(PYTHON_BIN)-config --includes 2>/dev/null)
PY_SUFFIX = $(shell $(PYTHON_BIN)-config --extension-suffix 2>/dev/null || $(PYTHON_BIN) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')")

CXXFLAGS = $(BASE_CXXFLAGS)

SOURCES = src/main.cpp src/beta_quadrature.cpp src/flamelet_io.cpp src/integration_core.cpp
CORE_SOURCES = src/beta_quadrature.cpp src/integration_core.cpp

.PHONY: all clean serial info check python

all: build/pdf_integrator python

build/pdf_integrator: $(SOURCES) | build
	$(CXX) $(CXXFLAGS) $(OMP_CFLAGS) $(INCLUDES) $(SOURCES) -o $@ $(LDFLAGS)

python: build/pdf_integrator_cpp$(PY_SUFFIX)

build/pdf_integrator_cpp$(PY_SUFFIX): src/bindings.cpp $(CORE_SOURCES) | build
	$(CXX) $(CXXFLAGS) $(OMP_CFLAGS) -shared -fPIC $(INCLUDES) $(PY_INCLUDES) $(PY_CFLAGS) src/bindings.cpp $(CORE_SOURCES) -o $@ -L$(OPENMP_LIB) $(OMP_LDFLAGS) $(PY_SHARED_LDFLAGS)

serial: build/pdf_integrator_serial

build/pdf_integrator_serial: $(SOURCES) | build
	$(CXX) $(CXXFLAGS) -DDISABLE_OPENMP $(INCLUDES) $(SOURCES) -o $@ -L$(HDF5_LIB) -lhdf5_cpp -lhdf5 -Wl,-rpath,$(HDF5_LIB)

build:
	mkdir -p build

clean:
	rm -rf build

info:
	@echo "PYTHON_BIN:    $(PYTHON_BIN)"
	@echo "CONDA_PREFIX:  $(CONDA_PREFIX)"
	@echo "HDF5 include:  $(HDF5_INC)"
	@echo "HDF5 lib:      $(HDF5_LIB)"
	@echo "OpenMP lib:    $(OPENMP_LIB)"

check:
	@echo "Checking dependencies under $(CONDA_PREFIX)..."
	@test -f $(HDF5_INC)/H5Cpp.h || (echo "ERROR: HDF5 C++ headers not found" && exit 1)
	@echo "HDF5 headers found."
	@if [ -f $(OPENMP_LIB)/$(OPENMP_LIB_FILE) ]; then \
		echo "OpenMP runtime found: $(OPENMP_LIB_FILE)"; \
	else \
		echo "WARNING: OpenMP runtime not found at $(OPENMP_LIB)/$(OPENMP_LIB_FILE)"; \
	fi
