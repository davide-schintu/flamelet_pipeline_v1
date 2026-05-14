#pragma once

#include <cmath>
#include <utility>
#include <vector>

namespace pdf {

/// Convert mean and variance to Beta distribution alpha/beta parameters
/// @param mean The mean of the distribution (0, 1)
/// @param variance The variance of the distribution
/// @return Pair of (alpha, beta) parameters
std::pair<double, double> mom2param(double mean, double variance);

/// Vectorized piecewise integration of properties over Beta PDF
/// Performs: integral[ Yi(xi) * BetaPDF(xi; alpha, beta) dxi ]
/// using piecewise linear interpolation of Yi
///
/// @param xi Z grid points (N points)
/// @param log_xi Precomputed log(Z) grid (N points)
/// @param log_1_xi Precomputed log(1-Z) grid (N points)
/// @param Yi Property values at each Z point (N x nProps)
/// @param alpha Beta distribution alpha parameter
/// @param beta_param Beta distribution beta parameter
/// @return Integrated values for each property (nProps)
std::vector<double> bct_vectorized(
    const std::vector<double>& xi,
    const std::vector<double>& log_xi,
    const std::vector<double>& log_1_xi,
    const std::vector<std::vector<double>>& Yi,
    double alpha,
    double beta_param
);

/// Process a single flamelet: integrate all properties over Z variance grid
/// @param Z The mixture fraction grid
/// @param log_Z Precomputed log(Z)
/// @param log_1_Z Precomputed log(1-Z)
/// @param Y_all Property matrix (nZ x nProps)
/// @param S_grid Normalized variance grid
/// @return Integrated properties (nZ x nS x nProps)
std::vector<std::vector<std::vector<double>>> process_single_flamelet(
    const std::vector<double>& Z,
    const std::vector<double>& log_Z,
    const std::vector<double>& log_1_Z,
    const std::vector<std::vector<double>>& Y_all,
    const std::vector<double>& S_grid
);

} // namespace pdf
