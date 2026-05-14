#include "../include/mixingRules.h"
#include <cmath>
#include <stdexcept>
// #include <Eigen/Dense>

void SRK_abc(double Tc, double Pc, double omega, double eps_kb, double sigma,
             double& ac, double& bc, double& cc)
{
    constexpr double k1 = 0.08664;
    constexpr double k2 = 0.42747;
    constexpr double Na = 6.02214076e23;
    constexpr double kb = 1.3806488e-23;

    if (omega == 0.0) {
        sigma *= 1e-10; // Converti Å in metri
        double eps = eps_kb * kb;
        cc = 0.48508;
        bc = 0.855 * Na * (sigma * sigma * sigma);
        ac = 5.55 * (Na * Na) * eps * (sigma * sigma * sigma);
    } else {
        cc = 0.48508 + 1.555171 * omega - 0.15613 * omega * omega;
        bc = (k1 * R0 * Tc) / Pc;
        ac = (k2 * R0 * R0 * Tc * Tc) / Pc;
    }
}

inline void alphaFunction(double T, double Tc, double cc, const std::string& EoS,
                   double& alpha, double& dalphadT, double& ddalphadT)
{
    // if (EoS == "PR" || EoS == "SRK") {
        double sqrt_T_Tc = std::sqrt(T / Tc);
        double term = 1.0 + cc * (1.0 - sqrt_T_Tc);
        alpha = term * term;
        dalphadT = -cc * std::sqrt(alpha) / std::sqrt(Tc * T);
        
        double T3_Tc_sqrt = std::sqrt(T * T * T * Tc);
        ddalphadT = (cc * cc) / (2.0 * T * Tc) + 
                    (cc / (2.0 * T3_Tc_sqrt)) * term;
    // }
    // else if (EoS == "RKPR") {
    //     double denom = 2.0 + T / Tc;
    //     alpha = std::pow(3.0 / denom, cc);
    //     dalphadT = -(cc * std::pow(3.0, cc)) / (Tc * std::pow(denom, cc + 1));
    //     ddalphadT = (std::pow(3.0, cc) * cc * (cc + 1)) / (Tc * Tc * std::pow(denom, cc + 2));
    // }
    // else {
    //     throw std::invalid_argument("The selected EoS is not available");
    // }
}

void VdWmixing_CR1(const std::vector<double>& X, double T,
                   const std::vector<Species>& Species, const std::string& EoS, double kij,
                   double& am, double& bm, double& dadTm, double& ddadTm)
{
    const size_t Ns = X.size();
    if (Species.size() != Ns)
        throw std::invalid_argument("Species vector and mole fractions size mismatch");

    // Preallocate and compute per-species constants
    std::vector<double> ac(Ns), bc(Ns), cc(Ns);
    std::vector<double> alpha(Ns), dalphadT(Ns), ddalphadT(Ns);
    std::vector<double> sqrt_ac(Ns), sqrt_alpha(Ns), inv_sqrt_alpha(Ns);

    for (size_t k = 0; k < Ns; ++k) {
        SRK_abc(Species[k].Tc, Species[k].Pc, Species[k].omega, Species[k].eps_kb, Species[k].sigma,
                ac[k], bc[k], cc[k]);
        alphaFunction(T, Species[k].Tc, cc[k], EoS, alpha[k], dalphadT[k], ddalphadT[k]);
        
        // Precompute expensive sqrt operations
        sqrt_ac[k] = std::sqrt(ac[k]);
        sqrt_alpha[k] = std::sqrt(alpha[k]);
        inv_sqrt_alpha[k] = 1.0 / sqrt_alpha[k];
    }

    // Covolume della miscela
    bm = 0.0;
    for (size_t i = 0; i < Ns; ++i)
        bm += X[i] * bc[i];

    // Initialize outputs
    am = 0.0;
    dadTm = 0.0;
    ddadTm = 0.0;

    const double one_minus_kij = 1.0 - kij;

    // Exploit symmetry: iterate j from i to Ns, apply factor 2 for off-diagonal
    for (size_t i = 0; i < Ns; ++i) {
        const double Xi = X[i];
        if (Xi < 1e-12) continue;

        const double sqrt_aci = sqrt_ac[i];
        const double sqrt_alphai = sqrt_alpha[i];
        const double inv_sqrt_alphai = inv_sqrt_alpha[i];
        const double dalphadTi = dalphadT[i];
        const double ddalphadT_i = ddalphadT[i];

        for (size_t j = i; j < Ns; ++j) {
            const double Xj = X[j];
            if (Xj < 1e-12) continue;

            const double sqrt_acj = sqrt_ac[j];
            const double sqrt_alphaj = sqrt_alpha[j];
            const double inv_sqrt_alphaj = inv_sqrt_alpha[j];
            const double dalphadTj = dalphadT[j];
            const double ddalphadT_j = ddalphadT[j];

            // Precomputed products
            const double sqrt_aci_acj = sqrt_aci * sqrt_acj;
            const double sqrt_alphai_alphaj = sqrt_alphai * sqrt_alphaj;
            const double inv_sqrt_alphai_alphaj = inv_sqrt_alphai * inv_sqrt_alphaj;

            const double term_x = Xi * Xj * sqrt_aci_acj;
            const double symmetry_factor = (i == j) ? 1.0 : 2.0;

            // am contribution
            am += symmetry_factor * term_x * sqrt_alphai_alphaj * one_minus_kij;

            // dadTm contribution (uses sqrt(alpha_i/alpha_j) = sqrt_alphai * inv_sqrt_alphaj)
            const double ratio_i_over_j = sqrt_alphai * inv_sqrt_alphaj;
            const double ratio_j_over_i = sqrt_alphaj * inv_sqrt_alphai;
            
            dadTm += symmetry_factor * term_x * (
                0.5 * ratio_i_over_j * dalphadTj +
                0.5 * ratio_j_over_i * dalphadTi
            );

            // ddadTm contribution
            const double inv_alphaj_32 = inv_sqrt_alphaj * inv_sqrt_alphaj * inv_sqrt_alphaj;
            const double inv_alphai_32 = inv_sqrt_alphai * inv_sqrt_alphai * inv_sqrt_alphai;
            
            ddadTm += symmetry_factor * term_x * (
                0.5 * inv_sqrt_alphai_alphaj * dalphadTi * dalphadTj
                - 0.25 * sqrt_alphai * inv_alphaj_32 * dalphadTj * dalphadTj
                - 0.25 * sqrt_alphaj * inv_alphai_32 * dalphadTi * dalphadTi
                + 0.5 * ratio_i_over_j * ddalphadT_j
                + 0.5 * ratio_j_over_i * ddalphadT_i
            );
        }
    }
}
void VdWmixing_CR1_ambm(const std::vector<double>& X, double T,
                        const std::vector<Species>& Species, const std::string& EoS, double kij,
                        double& am, double& bm)
{
    const size_t Ns = X.size();
    if (Species.size() != Ns)
        throw std::invalid_argument("Species vector and mole fractions size mismatch");

    std::vector<double> ac(Ns), bc(Ns), cc(Ns), alpha(Ns), sqrt_ac(Ns), sqrt_alpha(Ns);

    // if (EoS != "SRK")
    //     throw std::invalid_argument("The selected EoS is not available");

    // Precalcolo costanti
    for (size_t k = 0; k < Ns; ++k) {
        SRK_abc(Species[k].Tc, Species[k].Pc, Species[k].omega, Species[k].eps_kb, Species[k].sigma,
                ac[k], bc[k], cc[k]);
        double sqrt_T_Tc = std::sqrt(T / Species[k].Tc);
        alpha[k] = std::pow(1.0 + cc[k] * (1.0 - sqrt_T_Tc), 2);

        sqrt_ac[k] = std::sqrt(ac[k]);
        sqrt_alpha[k] = std::sqrt(alpha[k]);
    }

    // Calcolo bm
    bm = 0.0;
    for (size_t i = 0; i < Ns; ++i) {
        bm += X[i] * bc[i];
    }

    // Calcolo am con simmetria
    am = 0.0;
    for (size_t i = 0; i < Ns; ++i) {
        const double Xi = X[i];
        const double sqrt_aci_alphai = sqrt_ac[i] * sqrt_alpha[i];
        for (size_t j = i; j < Ns; ++j) {
            const double term = sqrt_aci_alphai * sqrt_ac[j] * sqrt_alpha[j];
            const double contribution = Xi * X[j] * term * (1.0 - kij);

            am += (i == j) ? contribution : 2.0 * contribution;
        }
    }
}

void VdWmixing_CR2(const std::vector<double>& X, double T,
                   const std::vector<Species>& Species, const std::string& EoS, double kij,
                   double& am, double& bm, double& dadTm, double& ddadTm)
{
    const size_t Ns = X.size();
    if (Species.size() != Ns)
        throw std::invalid_argument("Species vector and mole fractions size mismatch");

    std::vector<double> ac(Ns), bc(Ns), cc(Ns);
    for (size_t k = 0; k < Ns; ++k) {
        SRK_abc(Species[k].Tc, Species[k].Pc, Species[k].omega, Species[k].eps_kb, Species[k].sigma,
                ac[k], bc[k], cc[k]);
    }

    am = 0.0; 
    bm = 0.0;          
    dadTm = 0.0; 
    ddadTm = 0.0;

    for (size_t i = 0; i < Ns; ++i) {
        for (size_t j = 0; j < Ns; ++j) {
            double a_ij, alpha_ij, dalphadT_ij, ddalphadT_ij;
            if (i == j) {
                a_ij = ac[i];
                alphaFunction(T, Species[i].Tc, cc[i], EoS, alpha_ij, dalphadT_ij, ddalphadT_ij);
            } else {
                double vc_ij = 0.125 * std::pow(std::cbrt(Species[i].Vc) + std::cbrt(Species[j].Vc), 3);
                double Zc_ij = 0.5 * (Species[i].Zc + Species[j].Zc);
                double Tc_ij = std::sqrt(Species[i].Tc * Species[j].Tc) * (1.0 - kij);
                double pc_ij = Zc_ij * R0 * Tc_ij / vc_ij;
                double omega_ij = 0.5 * (Species[i].omega + Species[j].omega);

                double dummy_b, c_ij;
                SRK_abc(Tc_ij, pc_ij, omega_ij, 0.0, 0.0, a_ij, dummy_b, c_ij);
                alphaFunction(T, Tc_ij, c_ij, EoS, alpha_ij, dalphadT_ij, ddalphadT_ij);
            }
            double term = X[i] * X[j] * a_ij;
            am += term * alpha_ij;
            dadTm += term * dalphadT_ij;
            ddadTm += term * ddalphadT_ij;
        }
        bm += X[i] * bc[i]; 
    }
}

void VdWmixing_CR2_ambm(const std::vector<double>& X, double T, 
                   const std::vector<Species>& Species, const std::string& EoS, double kij, 
                   double& am, double& bm)
{
    const size_t Ns = X.size();
    std::vector<double> ac(Ns), bc(Ns), cc(Ns);
    for (size_t k = 0; k < Ns; ++k) {
        SRK_abc(Species[k].Tc, Species[k].Pc, Species[k].omega, Species[k].eps_kb, Species[k].sigma,
                ac[k], bc[k], cc[k]);
    }

    am = 0.0; bm = 0.0;
    for (size_t i = 0; i < Ns; ++i) {
        for (size_t j = 0; j < Ns; ++j) {
            double a_ij, alpha_ij, dummy1, dummy2;
            if (i == j) {
                a_ij = ac[i];
                alphaFunction(T, Species[i].Tc, cc[i], EoS, alpha_ij, dummy1, dummy2);
            } else {
                double vc_ij = 0.125 * std::pow(std::cbrt(Species[i].Vc) + std::cbrt(Species[j].Vc), 3);
                double Zc_ij = 0.5 * (Species[i].Zc + Species[j].Zc);
                double Tc_ij = std::sqrt(Species[i].Tc * Species[j].Tc) * (1.0 - kij);
                double pc_ij = Zc_ij * R0 * Tc_ij / vc_ij;
                double omega_ij = 0.5 * (Species[i].omega + Species[j].omega);

                double dummy_b, c_ij;
                SRK_abc(Tc_ij, pc_ij, omega_ij, 0.0, 0.0, a_ij, dummy_b, c_ij);
                alphaFunction(T, Tc_ij, c_ij, EoS, alpha_ij, dummy1, dummy2);
            }
            am += X[i] * X[j] * a_ij * alpha_ij;
        }
        bm += X[i] * bc[i];
    }
}

void VdWmixing_CR1_fugacity(const std::vector<double>& X, double T,
                            const std::vector<Species>& Species, const std::string& EoS, double kij,
                            double& am, double& bm, double& dadTm, double& ddadTm,
                            std::vector<double>& sum_ref_i, std::vector<double>& bc_out)
{
    const size_t Ns = X.size();
    std::vector<double> ac(Ns), bc(Ns), cc(Ns);
    std::vector<double> alpha(Ns), dalphadT(Ns), ddalphadT(Ns);

    for (size_t k = 0; k < Ns; ++k) {
        SRK_abc(Species[k].Tc, Species[k].Pc, Species[k].omega, Species[k].eps_kb, Species[k].sigma,
                ac[k], bc[k], cc[k]);
        alphaFunction(T, Species[k].Tc, cc[k], EoS, alpha[k], dalphadT[k], ddalphadT[k]);
    }

    bm = 0.0;
    for (size_t i = 0; i < Ns; ++i) bm += X[i] * bc[i];

    am = 0.0;
    dadTm = 0.0;
    ddadTm = 0.0;
    sum_ref_i.assign(Ns, 0.0);
    bc_out = bc;

    const double one_minus_kij = 1.0 - kij;

    for (size_t i = 0; i < Ns; ++i) {
        double ref_i = 0.0;
        for (size_t j = 0; j < Ns; ++j) {
            double a_ij = std::sqrt(ac[i] * ac[j]) * std::sqrt(alpha[i] * alpha[j]) * one_minus_kij;
            ref_i += X[j] * a_ij;

            double term_x = X[i] * X[j] * std::sqrt(ac[i] * ac[j]);
            am += term_x * std::sqrt(alpha[i] * alpha[j]) * one_minus_kij;

            dadTm += term_x * (0.5 * std::sqrt(alpha[i] / alpha[j]) * dalphadT[j] +
                               0.5 * std::sqrt(alpha[j] / alpha[i]) * dalphadT[i]) * one_minus_kij;

            ddadTm += term_x * ((1.0 / (2.0 * std::sqrt(alpha[i] * alpha[j]))) * dalphadT[i] * dalphadT[j]
                              - (0.25 * std::sqrt(alpha[i] / std::pow(alpha[j], 3))) * (dalphadT[j] * dalphadT[j])
                              - (0.25 * std::sqrt(alpha[j] / std::pow(alpha[i], 3))) * (dalphadT[i] * dalphadT[i])
                              + (0.5 * std::sqrt(alpha[i] / alpha[j])) * ddalphadT[j]
                              + (0.5 * std::sqrt(alpha[j] / alpha[i])) * ddalphadT[i]) * one_minus_kij;
        }
        sum_ref_i[i] = ref_i;
    }
}

void VdWmixing_CR2_fugacity(const std::vector<double>& X, double T,
                            const std::vector<Species>& Species, const std::string& EoS, double kij,
                            double& am, double& bm, double& dadTm, double& ddadTm,
                            std::vector<double>& sum_ref_i, std::vector<double>& bc_out)
{
    const size_t Ns = X.size();
    std::vector<double> ac(Ns), bc(Ns), cc(Ns);
    for (size_t k = 0; k < Ns; ++k) {
        SRK_abc(Species[k].Tc, Species[k].Pc, Species[k].omega, Species[k].eps_kb, Species[k].sigma,
                ac[k], bc[k], cc[k]);
    }

    am = 0.0; bm = 0.0; dadTm = 0.0; ddadTm = 0.0;
    sum_ref_i.assign(Ns, 0.0);
    bc_out = bc;

    for (size_t i = 0; i < Ns; ++i) {
        double ref_i = 0.0;
        for (size_t j = 0; j < Ns; ++j) {
            double a_ij, alpha_ij, dalphadT_ij, ddalphadT_ij;
            if (i == j) {
                a_ij = ac[i];
                alphaFunction(T, Species[i].Tc, cc[i], EoS, alpha_ij, dalphadT_ij, ddalphadT_ij);
            } else {
                double vc_ij = 0.125 * std::pow(std::cbrt(Species[i].Vc) + std::cbrt(Species[j].Vc), 3);
                double Zc_ij = 0.5 * (Species[i].Zc + Species[j].Zc);
                double Tc_ij = std::sqrt(Species[i].Tc * Species[j].Tc) * (1.0 - kij);
                double pc_ij = Zc_ij * R0 * Tc_ij / vc_ij;
                double omega_ij = 0.5 * (Species[i].omega + Species[j].omega);

                double dummy_b, c_ij;
                SRK_abc(Tc_ij, pc_ij, omega_ij, 0.0, 0.0, a_ij, dummy_b, c_ij);
                alphaFunction(T, Tc_ij, c_ij, EoS, alpha_ij, dalphadT_ij, ddalphadT_ij);
            }
            double term = X[i] * X[j] * a_ij;
            am += term * alpha_ij;
            dadTm += term * dalphadT_ij;
            ddadTm += term * ddalphadT_ij;
            ref_i += X[j] * a_ij * alpha_ij;
        }
        sum_ref_i[i] = ref_i;
        bm += X[i] * bc[i];
    }
}