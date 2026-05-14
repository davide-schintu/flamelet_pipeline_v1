#ifndef BETA_QUADRATURE_HPP
#define BETA_QUADRATURE_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <string>

namespace pdf {

class BetaQuadrature {
public:
    BetaQuadrature();
    
    // Set mean and variance to compute weights
    void Set(double mean, double variance);
    
    // Compute integral of interpolated function
    // Integral = Sum(f[i] * weight[i]) / BetaInf
    double IntegralNormalized(const std::vector<double>& f_b, 
                              const std::vector<double>& f_f, 
                              const std::vector<double>& f_center);
    
    // Accessors for grid points (needed for interpolation setup)
    const std::vector<double>& get_xb() const { return xb_; }
    const std::vector<double>& get_xf() const { return xf_; }
    const std::vector<double>& get_xcenter() const { return xcenter_; }
    
    size_t get_Nsub() const { return Nsub_; }
    size_t get_N() const { return N_; }
    size_t get_M() const { return M_; }

private:
    double FlatIntegral();
    void ErrorMessage(const std::string& message);

    // Configuration
    double epsilon_;
    unsigned int Nsub_;
    unsigned int N_;
    unsigned int M_; // Center points

    // Grids
    std::vector<double> dhb_, dhf_;
    std::vector<double> xb_, xf_, xcenter_;
    double dhcenter_;

    // PDF values (weights)
    std::vector<double> yb_, yf_, ycenter_;
    
    // Parameters
    double a_, b_;
    double extreme_a_, extreme_b_;
    double BetaInf_;
};

} // namespace pdf

#endif // BETA_QUADRATURE_HPP
