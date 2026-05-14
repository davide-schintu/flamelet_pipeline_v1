#include "integration_core.hpp"
#include "beta_quadrature.hpp"
#include "interpolation.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#ifndef DISABLE_OPENMP
#include <omp.h>
#endif

namespace core {

std::vector<double> integrate_property(
    const std::vector<double>& Z_grid,
    const std::vector<double>& S_grid,
    const std::vector<double>& flattened_input, // [nP * nChi * nZ]
    size_t nP,
    size_t nChi
) {
    size_t nZ = Z_grid.size();
    size_t nS = S_grid.size();
    
    // Output size: P * Chi * Z * S
    size_t total_output = nP * nChi * nZ * nS;
    std::vector<double> output(total_output, 0.0);
    
    // Initialize Master Quadrature and Interpolators just once
    pdf::BetaQuadrature master_quad; // Grids are fixed
    
    pdf::FastInterpolator interp_b, interp_f, interp_center;
    interp_b.Initialize(Z_grid, master_quad.get_xb());
    interp_f.Initialize(Z_grid, master_quad.get_xf());
    interp_center.Initialize(Z_grid, master_quad.get_xcenter());
    
#ifndef DISABLE_OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (size_t i = 0; i < nP; ++i) {
        for (size_t j = 0; j < nChi; ++j) {
            
            // Extract profile Y(Z) for this (P, Chi)
            size_t base_idx_3d = i * nChi * nZ + j * nZ;
            
            // Scratch arrays (thread-local)
            std::vector<double> y_prop(nZ);
            for (size_t k = 0; k < nZ; ++k) {
                y_prop[k] = flattened_input[base_idx_3d + k];
            }
            
            // Check max value to skip empty flamelets (optional optimization)
            double max_val = 0.0;
            for(double v : y_prop) max_val = std::max(max_val, v);
            // If max_val is very small, we might skip, but better be safe and compute
            
            // Local Quadrature object (cheap to construct, holds weights)
            pdf::BetaQuadrature local_quad; // Re-constructs grid (cheap) or we could clone?
            // Actually BetaQuadrature construction involves a few loops (generating grid).
            // It might be better to allow "setting" grid from master.
            // But 300 iterations is negligible compared to the integration work.
            // Let's optimize if needed. For now constructor is fine.
            
            // Interpolate ONCE per flamelet? 
            // NO. The weights depend on Variance. The interpolation depends on Y (flamelet).
            // Wait. OpenSMOKE interpolates Y onto the quadrature grid.
            // The quadrature grid is static. 
            // So we can interpolate Y -> Y_quad ONCE per flamelet!
            // Then for each S (variance), we just sum Y_quad * weights.
            
            // Optimization: Interpolate outside the S loop!
            // -------------------------------------------------------------
            std::vector<double> f_b(local_quad.get_xb().size());
            std::vector<double> f_f(local_quad.get_xf().size());
            std::vector<double> f_center(local_quad.get_xcenter().size());
            
            interp_b.Interpolate(y_prop, f_b);
            interp_f.Interpolate(y_prop, f_f);
            interp_center.Interpolate(y_prop, f_center);
            
            // Loop over Variance (S)
            for (size_t is = 0; is < nS; ++is) {
                double s = S_grid[is];
                // Loop over iz (Z_mean location)
                
                for (size_t iz = 0; iz < nZ; ++iz) {
                    double mean = Z_grid[iz];
                    double mazvard = mean * (1.0 - mean);
                    double variance = mazvard * s;
                    
                    // Laminar / Boundaries
                    size_t out_idx = i * nChi * nZ * nS 
                                   + j * nZ * nS 
                                   + iz * nS + is;

                    if (s <= 1e-6) {
                        output[out_idx] = y_prop[iz];
                        continue;
                    }
                    if (s >= 1.0 - 1e-6) {
                        output[out_idx] = y_prop[nZ-1]*mean + y_prop[0]*(1.0-mean);
                        continue;
                    }
                    if (mean < 1e-6) { output[out_idx] = y_prop[0]; continue; }
                    if (mean > 1.0-1e-6) { output[out_idx] = y_prop[nZ-1]; continue; }

                    // Set weights for this (mean, var) pair
                    local_quad.Set(mean, variance);
                    
                    // Integrate using PRE-INTERPOLATED values?
                    // NO! The interpolation maps Y(Z) -> Y(quad_grid).
                    // Does the quadrature grid depend on (mean, var)?
                    // OpenSMOKE: "Nsub_ = static_cast<unsigned int>(std::log10(0.1 / epsilon_));"
                    // The grid construction in constructor depends on epsilon_ (constant).
                    // The grid points xb_, xf_ are FIXED in [0,1].
                    // So YES, the grid is static. 
                    // So YES, we can interpolate Y once!
                    
                    // Wait, earlier I said "Interpolate Y -> Integrate".
                    // If Y is interpolated once, then we reuse f_b, f_f, f_center for all mean/var combinations?
                    // Yes! Because Y doesn't change. The WEIGHTS change with mean/var.
                    
                    // This is a HUGE optimization I missed in main.cpp!
                    // In main.cpp I was interpolating inside the p loop?
                    // No, inside the is_var loop I was interpolating?
                    // Let's check main.cpp logic.
                    // main.cpp:
                    // loop iz:
                    //   loop s:
                    //     Set(mean, var)
                    //     loop p:
                    //       interpolate(Y[p], f) -> integrate
                    
                    // In main.cpp, we perform interpolation for EACH (mean, var) pair?
                    // Can we hoist it out?
                    // The interpolation maps Z_grid -> Quad_grid. 
                    // Quad_grid is fixed.
                    // So yes, interpolation can be done once per flamelet property!
                    
                    // In this function, we process one property at a time.
                    // So for this flamelet (i,j) and property y_prop:
                    // 1. Interpolate y_prop -> f_b, f_f, f_center (ONCE)
                    // 2. Loop iz, Loop is:
                    //      Set weights (depends on mean, var)
                    //      Sum product
                    
                    // Verify main.cpp logic:
                    /*
                    for (size_t iz = 0; iz < nZ; ++iz) {
                       ...
                       for (size_t is_var = 0; is_var < nS; ++is_var) {
                          ...
                          quad.Set(mean, variance);
                          for (size_t p = 0; p < nProps; ++p) {
                             interp...
                             integrate...
                          }
                       }
                    }
                    */
                    // In main.cpp, interpolation was inside the innermost loop!
                    // This was inefficient! 
                    // But wait, `interp.Interpolate` uses `indices` and `weights` which are precomputed.
                    // So it's just a multiply-add loop.
                    // Still, doing it N_Z * N_S times is redundant.
                    // It should be done ONCE per flamelet.
                    
                    // I will implement this OPTIMIZATION here.
                    
                    output[out_idx] = local_quad.IntegralNormalized(f_b, f_f, f_center);
                }
            }
        }
    }
    
    return output;
}

std::vector<double> integrate_density_opensmoke(
    const std::vector<double>& Z_grid,
    const std::vector<double>& S_grid,
    const std::vector<double>& rho,
    size_t nP,
    size_t nChi
) {
    std::vector<double> inv_rho(rho.size(), 0.0);
    for (size_t i = 0; i < rho.size(); ++i) {
        if (rho[i] <= 0.0) {
            throw std::runtime_error("rho must be positive for OpenSMOKE density convolution");
        }
        inv_rho[i] = 1.0 / rho[i];
    }

    std::vector<double> inv_rho_mean = integrate_property(Z_grid, S_grid, inv_rho, nP, nChi);
    std::vector<double> rho_mean(inv_rho_mean.size(), 0.0);
    for (size_t i = 0; i < inv_rho_mean.size(); ++i) {
        if (inv_rho_mean[i] <= 0.0) {
            throw std::runtime_error("convoluted inverse density must be positive");
        }
        rho_mean[i] = 1.0 / inv_rho_mean[i];
    }
    return rho_mean;
}

std::vector<double> integrate_favre_property(
    const std::vector<double>& Z_grid,
    const std::vector<double>& S_grid,
    const std::vector<double>& rho,
    const std::vector<double>& property,
    const std::vector<double>& rho_mean,
    size_t nP,
    size_t nChi
) {
    if (rho.size() != property.size()) {
        throw std::runtime_error("rho and property size mismatch");
    }

    std::vector<double> rho_property(property.size(), 0.0);
    for (size_t i = 0; i < property.size(); ++i) {
        rho_property[i] = rho[i] * property[i];
    }

    std::vector<double> favre = integrate_property(Z_grid, S_grid, rho_property, nP, nChi);
    if (favre.size() != rho_mean.size()) {
        throw std::runtime_error("rho_mean and integrated property size mismatch");
    }

    for (size_t i = 0; i < favre.size(); ++i) {
        if (rho_mean[i] <= 0.0) {
            throw std::runtime_error("rho_mean must be positive for Favre convolution");
        }
        favre[i] /= rho_mean[i];
    }
    return favre;
}

} // namespace core
