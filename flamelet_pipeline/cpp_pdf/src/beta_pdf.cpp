#include "beta_pdf.hpp"
#include <boost/math/special_functions/beta.hpp>
#include <algorithm>
#include <cmath>

namespace pdf {

std::pair<double, double> mom2param(double m, double v) {
    constexpr double epsilon = 1e-10;
    constexpr double tiny = 1e-32;
    
    // Boundary / Invalid check - same logic as Python
    if (m <= epsilon || m >= 1.0 - epsilon || 
        v <= epsilon || v >= (m * (1.0 - m)) - epsilon) {
        return {tiny, tiny};
    }
    
    // Calculate alpha (a) and beta (b)
    // v = m(1-m)/(a+b+1) -> a+b = m(1-m)/v - 1
    double term = (m * (1.0 - m) / v) - 1.0;
    if (term < epsilon) {
        term = epsilon;
    }
    
    double a = m * term;
    double b = (1.0 - m) * term;
    
    return {a, b};
}

std::vector<double> bct_vectorized(
    const std::vector<double>& xi,
    const std::vector<double>& log_xi,
    const std::vector<double>& log_1_xi,
    const std::vector<std::vector<double>>& Yi,
    double alpha,
    double beta_param
) {
    const size_t N = xi.size();
    const size_t nProps = Yi.empty() ? 0 : Yi[0].size();
    
    std::vector<double> result(nProps, 0.0);
    
    if (N < 2 || nProps == 0) {
        return result;
    }

    // Precompute beta factors
    double inv_beta = 1.0 / boost::math::beta(alpha, beta_param);
    double beta_fact_1 = boost::math::beta(alpha + 1.0, beta_param) * inv_beta;
    
    // Term for recurrence: I_x(a+1, b) = I_x(a, b) - x^a(1-x)^b / (a*B(a,b))
    // x^a(1-x)^b can be computed as exp(a*log(x) + b*log(1-x))
    double recurrence_factor = inv_beta / alpha;

    // Process each segment
    for (size_t seg = 0; seg < N - 1; ++seg) {
        double x1 = xi[seg];
        double x2 = xi[seg + 1];
        double dx = x2 - x1;
        
        // Filter small dx to avoid division by zero
        if (std::abs(dx) < 1e-15) {
            continue;
        }
        
        // Beta function integrals
        // Use recurrence relation to save calls
        
        double I0_x1 = boost::math::ibeta(alpha, beta_param, x1);
        double I1_x1;
        
        if (x1 > 1e-9 && x1 < 1.0 - 1e-9) {
            // Use precomputed logs: exp(alpha * log(x) + beta * log(1-x))
            double log_term = alpha * log_xi[seg] + beta_param * log_1_xi[seg];
            double term = std::exp(log_term) * recurrence_factor;
            I1_x1 = I0_x1 - term;
        } else {
            // Near boundaries, fall back to direct computation
            I1_x1 = boost::math::ibeta(alpha + 1.0, beta_param, x1);
        }

        double I0_x2 = boost::math::ibeta(alpha, beta_param, x2);
        double I1_x2;
        
        if (x2 > 1e-9 && x2 < 1.0 - 1e-9) {
            double log_term = alpha * log_xi[seg + 1] + beta_param * log_1_xi[seg + 1];
            double term = std::exp(log_term) * recurrence_factor;
            I1_x2 = I0_x2 - term;
        } else {
            I1_x2 = boost::math::ibeta(alpha + 1.0, beta_param, x2);
        }

        double Diff0 = I0_x2 - I0_x1;
        double Diff1 = I1_x2 - I1_x1;
        
        double Term_M = Diff1 * beta_fact_1;
        double Term_B = Diff0;
        
        // Accumulate contribution for each property
        for (size_t p = 0; p < nProps; ++p) {
            double y1 = Yi[seg][p];
            double y2 = Yi[seg + 1][p];
            
            // Linear segment coefficients (y = m*x + b)
            double m_coef = (y2 - y1) / dx;
            double b_int = y1 - m_coef * x1;
            
            result[p] += m_coef * Term_M + b_int * Term_B;
        }
    }
    
    return result;
}

std::vector<std::vector<std::vector<double>>> process_single_flamelet(
    const std::vector<double>& Z,
    const std::vector<double>& log_Z,
    const std::vector<double>& log_1_Z,
    const std::vector<std::vector<double>>& Y_all,
    const std::vector<double>& S_grid
) {
    const size_t nZ = Z.size();
    const size_t nS = S_grid.size();
    const size_t nProps = Y_all.empty() ? 0 : Y_all[0].size();
    
    // Output: [nZ][nS][nProps]
    std::vector<std::vector<std::vector<double>>> Out_all(
        nZ, std::vector<std::vector<double>>(
            nS, std::vector<double>(nProps, 0.0)
        )
    );
    
    // Iterate over target Zmean and variance S
    for (size_t iz = 0; iz < nZ; ++iz) {
        double z_mean = Z[iz];
        double mazvard = z_mean * (1.0 - z_mean);  // Maximum variance at this Z
        
        for (size_t is_var = 0; is_var < nS; ++is_var) {
            double s_var = S_grid[is_var];
            
            // Case 1: Laminar (S ~ 0)
            if (s_var <= 1e-6) {
                Out_all[iz][is_var] = Y_all[iz];
                continue;
            }
            
            // Case 2: Max Mixing (S ~ 1) -> Double Delta
            if (s_var >= 1.0 - 1e-6) {
                // Linear combination of pure oxidizer and fuel
                for (size_t p = 0; p < nProps; ++p) {
                    double val0 = Y_all[0][p];          // Pure Oxidizer (Z=0)
                    double val1 = Y_all[nZ - 1][p];     // Pure Fuel (Z=1)
                    Out_all[iz][is_var][p] = val1 * z_mean + val0 * (1.0 - z_mean);
                }
                continue;
            }
            
            // Case 3: Boundaries of Zmean space (0 or 1)
            if (z_mean < 1e-6) {
                Out_all[iz][is_var] = Y_all[0];
                continue;
            }
            if (z_mean > 1.0 - 1e-6) {
                Out_all[iz][is_var] = Y_all[nZ - 1];
                continue;
            }
            
            // Standard PDF Convolution
            double variance = mazvard * s_var;
            auto [a, b] = mom2param(z_mean, variance);
            
            if (a < 1e-6 || b < 1e-6) {
                // Fallback to laminar if params invalid
                Out_all[iz][is_var] = Y_all[iz];
            } else {
                // Perform Integration
                Out_all[iz][is_var] = bct_vectorized(Z, log_Z, log_1_Z, Y_all, a, b);
            }
        }
    }
    
    return Out_all;
}

} // namespace pdf
