#ifndef INTEGRATION_CORE_HPP
#define INTEGRATION_CORE_HPP

#include <vector>
#include <string>

namespace core {

// Core function to integrate a full dataset of flamelets
// Z_grid: [nZ]
// Y_data: [nP * nChi * nZ * nProperties] (Flattened or structured external access? Let's use flattened simple vector for easy binding)
// Actually, let's keep it structured:
// input_data: Vector of properties. Each property is [nP * nChi * nZ]
// 
// Returns: output_data [nP * nChi * nZ * nS * nProperties] (or similar)
//
// To keep it simple and efficient for binding:
// We'll perform integration on a single "property block" like T or alpha.
// integrate_property(Z_grid, S_grid, flattened_3d_property) -> flattened_4d_result

std::vector<double> integrate_property(
    const std::vector<double>& Z_grid,
    const std::vector<double>& S_grid,
    const std::vector<double>& flattened_input, // [nP * nChi * nZ]
    size_t nP,
    size_t nChi
);

std::vector<double> integrate_density_opensmoke(
    const std::vector<double>& Z_grid,
    const std::vector<double>& S_grid,
    const std::vector<double>& rho,
    size_t nP,
    size_t nChi
);

std::vector<double> integrate_favre_property(
    const std::vector<double>& Z_grid,
    const std::vector<double>& S_grid,
    const std::vector<double>& rho,
    const std::vector<double>& property,
    const std::vector<double>& rho_mean,
    size_t nP,
    size_t nChi
);

} // namespace core

#endif // INTEGRATION_CORE_HPP
