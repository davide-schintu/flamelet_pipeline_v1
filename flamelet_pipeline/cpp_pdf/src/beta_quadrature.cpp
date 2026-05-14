#include "beta_quadrature.hpp"
#include <algorithm>
#include <cmath>

namespace pdf {

BetaQuadrature::BetaQuadrature() {
    epsilon_ = 1.e-6;
    // Logarithmic sub-intervals count: log10(0.1 / 1e-6) = 5
    Nsub_ = static_cast<unsigned int>(std::log10(0.1 / epsilon_));
    N_ = 40;  // Increased from 20 for accuracy
    M_ = 200; // Increased from 100 for accuracy

    dhb_.resize(Nsub_, 0.0);
    dhf_.resize(Nsub_, 0.0);

    xb_.resize(Nsub_ * N_, 0.0);
    yb_.resize(Nsub_ * N_, 0.0);
    xf_.resize(Nsub_ * N_, 0.0);
    yf_.resize(Nsub_ * N_, 0.0);

    xcenter_.resize(M_, 0.0);
    ycenter_.resize(M_, 0.0);

    // 1. Center interval [0.1, 0.9]
    dhcenter_ = (0.9 - 0.1) / double(M_);
    xcenter_[0] = 0.1 + dhcenter_ * 0.50;
    for (unsigned int j = 1; j < M_; j++)
        xcenter_[j] = xcenter_[j - 1] + dhcenter_;

    // 2. Backward Interval [0, 0.1] (Logarithmic clustering near 0)
    for (unsigned int k = 0; k < Nsub_; k++)
        dhb_[k] = epsilon_ * (std::pow(10.0, double(k + 1)) - std::pow(10.0, double(k))) / double(N_);

    for (unsigned int k = 0; k < Nsub_; k++)
        xb_[k * N_] = epsilon_ * std::pow(10.0, double(k)) + dhb_[k] * 0.50;

    for (unsigned int k = 0; k < Nsub_; k++)
        for (unsigned int j = k * N_ + 1; j < (k + 1) * N_; j++)
            xb_[j] = xb_[j - 1] + dhb_[k];

    // 3. Forward Interval [0.9, 1.0] (Symmetric to Backward)
    dhf_ = dhb_;
    for (unsigned int j = 0; j < Nsub_ * N_; j++)
        xf_[j] = 1.0 - xb_[j];
}

void BetaQuadrature::ErrorMessage(const std::string& message) {
    std::cerr << "\nClass: BetaQuadrature\nError: " << message << "\n";
    std::exit(-1);
}

void BetaQuadrature::Set(double mean, double variance) {
    // Calculate alpha (a) and beta (b) parameters
    double tmp = mean * (1.0 - mean) / variance - 1.0;
    a_ = mean * tmp;
    b_ = (1.0 - mean) * tmp;

    if (a_ <= 0.0) ErrorMessage("The a coefficient must be positive");
    if (b_ <= 0.0) ErrorMessage("The b coefficient must be positive");

    // Clipping for very large coefficients (to avoid overflow/instability)
    if (a_ > 1.0 && b_ > 1.0) {
        const double fmax = 1.0 / (1.0 + (b_ - 1.0) / (a_ - 1.0));

        if (a_ > 500.0) {
            a_ = 500.0;
            b_ = (a_ - 1.0 - fmax * (a_ - 2.0)) / fmax;
        } else if (b_ > 500.0) {
            b_ = 500.0;
            a_ = (1.0 + fmax * (b_ - 2.0)) / (1.0 - fmax);
        }
    }

    const double a_minus_one = a_ - 1.0;
    const double b_minus_one = b_ - 1.0;

    // Compute PDF values at quadrature points
    for (unsigned int j = 0; j < Nsub_ * N_; j++)
        yb_[j] = std::pow(xb_[j], a_minus_one) * std::pow(1.0 - xb_[j], b_minus_one); // Note: 1-xb is computed directly? No, xb is small.
        // Wait, OpenSMOKE used: yb_(j) = pow(xb, a-1) * pow(xf, b-1)? 
        // In OpenSMOKE xf[j] = 1 - xb[j]. But xb and xf arrays are separate.
        // Let's stick to explicit math: pow(x, a-1)*pow(1-x, b-1).
        // Since xb is near 0, 1-xb is near 1.
    
    // Correction: In OpenSMOKE code:
    // yb_(j) = std::pow(xb_(j), a_minus_one)*std::pow(xf_(j), b_minus_one); 
    // This is WRONG if xf_ is the forward grid. It meant (1-xb).
    // Ah, wait. OpenSMOKE defines xf[j] as 1 - xb[j].
    // BUT xf_ array stores points near 1.
    // So for yb calculation, we need (1-xb). Is xf_[j] equal to 1-xb_[j]?
    // Yes: xf_[j] = 1. - xb_[j] in loop line 86.
    // So yes, using xf_[j] is correct for (1-xb).
    
    for (unsigned int j = 0; j < Nsub_ * N_; j++)
        yb_[j] = std::pow(xb_[j], a_minus_one) * std::pow(xf_[j], b_minus_one);

    for (unsigned int j = 0; j < Nsub_ * N_; j++)
        yf_[j] = std::pow(xf_[j], a_minus_one) * std::pow(xb_[j], b_minus_one);

    for (unsigned int j = 0; j < M_; j++)
        ycenter_[j] = std::pow(xcenter_[j], a_minus_one) * std::pow(1.0 - xcenter_[j], b_minus_one);

    // Contributions from singularity tails (analytic approximation)
    extreme_a_ = std::pow(epsilon_, a_) / a_;
    extreme_b_ = std::pow(epsilon_, b_) / b_;

    BetaInf_ = FlatIntegral();
}

double BetaQuadrature::FlatIntegral() {
    double sum = 0.0;

    // Backward grid integration
    for (unsigned int k = 0; k < Nsub_; k++) {
        double sum_partial = 0.0;
        for (unsigned int j = k * N_; j < (k + 1) * N_; j++)
            sum_partial += yb_[j];
        sum_partial *= dhb_[k];
        sum += sum_partial;
    }

    // Forward grid integration
    for (unsigned int k = 0; k < Nsub_; k++) {
        double sum_partial = 0.0;
        for (unsigned int j = k * N_; j < (k + 1) * N_; j++)
            sum_partial += yf_[j];
        sum_partial *= dhf_[k];
        sum += sum_partial;
    }

    // Center grid integration
    {
        double sum_partial = 0.0;
        for (unsigned int j = 0; j < M_; j++)
            sum_partial += ycenter_[j];
        sum_partial *= dhcenter_;
        sum += sum_partial;
    }

    sum += extreme_a_ + extreme_b_;

    return sum;
}

double BetaQuadrature::IntegralNormalized(const std::vector<double>& f_b, 
                                          const std::vector<double>& f_f, 
                                          const std::vector<double>& f_center) {
    double sum = 0.0;

    // Backward
    for (unsigned int k = 0; k < Nsub_; k++) {
        double sum_partial = 0.0;
        for (unsigned int j = k * N_; j < (k + 1) * N_; j++)
            sum_partial += yb_[j] * f_b[j];
        sum_partial *= dhb_[k];
        sum += sum_partial;
    }

    // Forward
    for (unsigned int k = 0; k < Nsub_; k++) {
        double sum_partial = 0.0;
        for (unsigned int j = k * N_; j < (k + 1) * N_; j++)
            sum_partial += yf_[j] * f_f[j];
        sum_partial *= dhf_[k];
        sum += sum_partial;
    }

    // Center
    {
        double sum_partial = 0.0;
        for (unsigned int j = 0; j < M_; j++)
            sum_partial += ycenter_[j] * f_center[j];
        sum_partial *= dhcenter_;
        sum += sum_partial;
    }
    
    // Tails: assume f is constant at boundaries f(0) and f(1)
    // f_b[0] is close to 0? epsilon/2.
    // Note: OpenSMOKE uses f_a_ and f_b_ which are f(0) and f(1).
    // How do we get them? We assume the caller provides correct vectors or we need f(0)/f(1).
    // Let's assume f at singularity is f_b[0] (closest to 0) and f_f[0] (closest to 1).
    // Better: pass f(0) and f(1) separately? 
    // Or just use the closest point. xb[0] is very small (1e-6/20?).
    // OpenSMOKE passes interpolated f(0) and f(1) explicitly.
    // I will use f_b[0] approx f(0) and f_f[0] approx f(1) for simplicity, 
    // or better, I should ask for f_start/f_end.
    // Actually, f_b[0] is at x=0.5*epsilon/N.
    // Let's use f_b[0] and f_f[0] as reasonable approximations.
    
    double val_0 = f_b.empty() ? 0.0 : f_b[0]; 
    double val_1 = f_f.empty() ? 0.0 : f_f[0]; 

    sum += extreme_a_ * val_0 + extreme_b_ * val_1;

    return sum / BetaInf_;
}

} // namespace pdf
