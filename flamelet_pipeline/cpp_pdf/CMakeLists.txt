cmake_minimum_required(VERSION 3.16)
project(pdf_integrator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PDF_ENABLE_OPENMP "Enable OpenMP parallelism" ON)
option(PDF_ENABLE_NATIVE_OPT "Enable -march=native optimizations" OFF)
option(PDF_ENABLE_FAST_MATH "Enable -ffast-math optimizations" OFF)

find_package(HDF5 REQUIRED COMPONENTS CXX)

set(PDF_COMMON_SOURCES
    src/main.cpp
    src/flamelet_io.cpp
    src/integration_core.cpp
    src/beta_quadrature.cpp
)

set(PDF_COMMON_INCLUDES
    include
    ${HDF5_INCLUDE_DIRS}
)

set(PDF_RELEASE_FLAGS -O3)
if(PDF_ENABLE_NATIVE_OPT)
    list(APPEND PDF_RELEASE_FLAGS -march=native)
endif()
if(PDF_ENABLE_FAST_MATH)
    list(APPEND PDF_RELEASE_FLAGS -ffast-math)
endif()

add_executable(pdf_integrator ${PDF_COMMON_SOURCES})
target_include_directories(pdf_integrator PRIVATE ${PDF_COMMON_INCLUDES})
target_link_libraries(pdf_integrator PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_compile_options(pdf_integrator PRIVATE
    $<$<CONFIG:Release>:${PDF_RELEASE_FLAGS}>
)

if(PDF_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(pdf_integrator PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(STATUS "OpenMP not found: building without OpenMP support")
        target_compile_definitions(pdf_integrator PRIVATE DISABLE_OPENMP)
    endif()
else()
    target_compile_definitions(pdf_integrator PRIVATE DISABLE_OPENMP)
endif()

add_executable(pdf_integrator_serial ${PDF_COMMON_SOURCES})
target_include_directories(pdf_integrator_serial PRIVATE ${PDF_COMMON_INCLUDES})
target_link_libraries(pdf_integrator_serial PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_compile_definitions(pdf_integrator_serial PRIVATE DISABLE_OPENMP)
target_compile_options(pdf_integrator_serial PRIVATE
    $<$<CONFIG:Release>:${PDF_RELEASE_FLAGS}>
)
