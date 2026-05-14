#pragma once

#include <string>
#include <vector>
#include <H5Cpp.h>

namespace io {

/// Structure holding flamelet database dimensions and data
struct FlameletData {
    // Grid coordinates
    std::vector<double> Z_grid;      // Mixture fraction (nZ)
    std::vector<double> S_grid;      // Normalized variance (nS) - output only
    std::vector<double> pressure;    // Pressure levels (nP)
    std::vector<double> chi;         // Scalar dissipation (nChi)
    
    // Dimensions
    size_t nP = 0;
    size_t nChi = 0;
    size_t nZ = 0;
    size_t nS = 0;
    
    // Property names to process
    std::vector<std::string> property_names;
    
    // 3D input data stored as [P][Chi][Z] flattened to 1D
    // Access: data[prop_idx][i * nChi * nZ + j * nZ + k]
    std::vector<std::vector<double>> input_data;
    
    // 4D output data stored as [P][Chi][Z][S] flattened
    std::vector<std::vector<double>> output_data;
};

/// Read input HDF5 flamelet database
/// @param filename Path to input .h5 file
/// @param property_names Properties to extract and process
/// @return FlameletData structure with loaded data
FlameletData readInputH5(
    const std::string& filename,
    const std::vector<std::string>& property_names
);

/// Write output HDF5 with integrated properties
/// @param filename Path to output .h5 file
/// @param data FlameletData with computed output_data
void writeOutputH5(
    const std::string& filename,
    const FlameletData& data
);

/// Copy base structure from input to output (Z, pressure, chi grids)
void copyBaseStructure(
    const std::string& input_file,
    const std::string& output_file,
    const std::vector<double>& S_grid
);

} // namespace io
