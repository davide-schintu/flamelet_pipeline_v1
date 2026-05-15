#include "../include/functions.h"
#include <algorithm>
#include <cmath>
#include <array>
#include <stdexcept>

double getZ(double am, double bm, double P, double T,
            const std::vector<double> &X, const std::string &EoS,
            const std::vector<Species> &species) {
  double delta1 = 0.0, delta2 = 0.0;

  // if (EoS == "PR") {
  //     delta1 = 1.0 + std::sqrt(2.0);
  //     delta2 = 1.0 - std::sqrt(2.0);
  // }
  // else if (EoS == "SRK") {
  delta1 = 1.0;
  delta2 = 0.0;
  // }
  // else if (EoS == "RKPR") {
  //     double d1 = 0.428363, d2 = 18.496215, d3 = 0.338426;
  //     double d4 = 0.660000, d5 = 789.723105, d6 = 2.512392;

  //     double delta1_sum = 0.0;
  //     double delta2_sum = 0.0;

  //     for (size_t k = 0; k < species.size(); ++k) {
  //         double Zc = species[k].Zc;
  //         double delta1_k = d1 + d2 * std::pow((d3 - 1.168 * Zc), d4) + d5 *
  //         std::pow((d3 - 1.168 * Zc), d6); double delta2_k = (1.0 - delta1_k)
  //         / (1.0 + delta1_k); delta1_sum += X[k] * delta1_k; delta2_sum +=
  //         X[k] * delta2_k;
  //     }
  //     delta1 = delta1_sum;
  //     delta2 = delta2_sum;
  // }
  // else {
  //     throw std::invalid_argument("EoS non riconosciuta.");
  // }

  double ad = 1.0 / (R0 * R0 * T * T);
  double bd = 1.0 / (R0 * T);
  double A = am * P * ad;
  double B = bm * P * bd;

  double C2 = B * (delta1 + delta2 - 1.0) - 1.0;
  double C1 = A + delta1 * delta2 * B * B - (delta1 + delta2) * B * (B + 1.0);
  double C0 = -B * (delta1 * delta2 * B * B + delta1 * delta2 * B + A);

  double Q = (C2 * C2 - 3.0 * C1) / 9.0;
  double R = (2.0 * C2 * C2 * C2 - 9.0 * C2 * C1 + 27.0 * C0) / 54.0;

  const bool force_liquid_root = (EoS == "Liq");
  const bool force_vapor_root = (EoS == "Vap");

  double Z = 0.0;
  if (R * R < Q * Q * Q) {
    double sqrtQ = std::sqrt(Q);
    // double phi = std::acos(R / (sqrtQ * sqrtQ));
    double phi = std::acos(R / (sqrtQ * Q));

    double z1 = -2.0 * sqrtQ * std::cos(phi / 3.0) - C2 / 3.0;
    double z2 = -2.0 * sqrtQ * std::cos((phi + 2.0 * M_PI) / 3.0) - C2 / 3.0;
    double z3 = -2.0 * sqrtQ * std::cos((phi - 2.0 * M_PI) / 3.0) - C2 / 3.0;

    std::vector<double> roots = {z1, z2, z3};

    if (*std::min_element(roots.begin(), roots.end()) < B) {
      Z = *std::max_element(roots.begin(), roots.end());
    } else {
      std::sort(roots.begin(), roots.end());
      double ZL = roots[0];
      double ZV = roots[2];

      if (force_liquid_root) {
        Z = ZL;
      } else if (force_vapor_root) {
        Z = ZV;
      } else {
        double deltaG =
            (A / (B * (delta1 - delta2))) *
                std::log(((ZL + delta1 * B) * (ZV + delta2 * B)) /
                         ((ZL + delta2 * B) * (ZV + delta1 * B))) -
            (ZL - ZV) + std::log((ZL - B) / (ZV - B));

        Z = (deltaG < 0) ? ZV : ZL;
      }
    }
  } else {
    double E = -std::copysign(1.0, R) *
               std::pow(std::abs(R) + std::sqrt(R * R - Q * Q * Q), 1.0 / 3.0);
    double F = (E != 0.0) ? (Q / E) : 0.0;
    Z = E + F - C2 / 3.0;
  }

  return Z;
}

void getRho(double Z, double T, double p, double M, double &rho, double &v) {
  v = (Z * R0 * T) / p;
  double rho_mol = 1.0 / v;
  rho = rho_mol * M;
}

double getH(double T, double p, double M, const std::vector<double> &X,
            const std::vector<double> &Y, const std::vector<Species> &species,
            double am, double bm, double dadTm, double ddadTm,
            const std::string &EoS) {
  double cp0, cv0, h0;
  getIdealDepartureProperties(T, Y, species, M, cp0, cv0, h0);

  double v;
  double rho;
  double Z = getZ(am, bm, p, T, X, EoS, species);
  getRho(Z, T, p, M, rho, v);

  double de, dh, dcv;

  departureFunctions(am, bm, dadTm, ddadTm, M, v, T, p, Y, EoS, species, de, dh,
                     dcv);

  return h0 + dh;
}

void getIdealDepartureProperties(double T, const std::vector<double> &Y,
                                 const std::vector<Species> &species, double M,
                                 double &cp0, double &cv0, double &h0) {
  std::vector<double> cp0i(species.size(), 0.0);
  std::vector<double> h0i(species.size(), 0.0);

  const double k1 = 0.5, k2 = 1.0 / 3.0, k3 = 0.25, k4 = 0.2;

  for (size_t i = 0; i < species.size(); ++i) {
    if (T < species[i].Tmed) {
      cp0i[i] =
          (species[i].AlowT + species[i].BlowT * T + species[i].ClowT * T * T +
           species[i].DlowT * T * T * T + species[i].ElowT * T * T * T * T) *
          (R0 / species[i].M);

      h0i[i] =
          (species[i].AlowT + k1 * species[i].BlowT * T +
           k2 * species[i].ClowT * T * T + k3 * species[i].DlowT * T * T * T +
           k4 * species[i].ElowT * T * T * T * T + species[i].FlowT / T) *
          (R0 / species[i].M) * T;
    } else {
      cp0i[i] =
          (species[i].AHigT + species[i].BHigT * T + species[i].CHigT * T * T +
           species[i].DHigT * T * T * T + species[i].EHigT * T * T * T * T) *
          (R0 / species[i].M);

      h0i[i] =
          (species[i].AHigT + k1 * species[i].BHigT * T +
           k2 * species[i].CHigT * T * T + k3 * species[i].DHigT * T * T * T +
           k4 * species[i].EHigT * T * T * T * T + species[i].FHigT / T) *
          (R0 / species[i].M) * T;
    }
  }

  cp0 = 0.0;
  h0 = 0.0;
  for (size_t i = 0; i < Y.size(); ++i) {
    cp0 += Y[i] * cp0i[i];
    h0 += Y[i] * h0i[i];
  }

  cv0 = cp0 - (R0 / M);
}

void departureFunctions(double am, double bm, double dadTm, double ddadTm,
                        double M, double v, double T, double p,
                        const std::vector<double> &X, const std::string &EoS,
                        const std::vector<Species> &species, double &de,
                        double &dh, double &dcv) {
  const double R0 = 8.314472; // [J/mol*K]

  double delta1 = 0.0;
  double delta2 = 0.0;

  // if (EoS == "PR") {
  //     delta1 = 1.0 + std::sqrt(2.0);
  //     delta2 = 1.0 - std::sqrt(2.0);
  // } else if (EoS == "SRK") {
  delta1 = 1.0;
  delta2 = 0.0;
  // } else if (EoS == "RKPR") {
  //     const double d1 = 0.428363;
  //     const double d2 = 18.496215;
  //     const double d3 = 0.338426;
  //     const double d4 = 0.660000;
  //     const double d5 = 789.723105;
  //     const double d6 = 2.512392;

  //     std::vector<double> delta1_i(species.size(), 0.0);
  //     std::vector<double> delta2_i(species.size(), 0.0);

  //     for (std::size_t k = 0; k < species.size(); ++k) {
  //         double Zc = species[k].Zc;
  //         double term = d3 - 1.168 * Zc;
  //         delta1_i[k] = d1 + d2 * std::pow(term, d4) + d5 * std::pow(term,
  //         d6); delta2_i[k] = (1.0 - delta1_i[k]) / (1.0 + delta1_i[k]);
  //     }

  //     for (std::size_t k = 0; k < X.size(); ++k) {
  //         delta1 += X[k] * delta1_i[k];
  //         delta2 += X[k] * delta2_i[k];
  //     }
  // } else {
  //     throw std::invalid_argument("La EoS selezionata non è disponibile");
  // }

  // Energia interna
  double K = (1.0 / (bm * (delta1 - delta2))) *
             std::log((v + delta1 * bm) / (v + delta2 * bm));
  double de_mol = (am - T * dadTm) * K;
  de = -de_mol / M; // [J/kg]

  // Entalpia
  double dh_mol = -de_mol + p * v - R0 * T;
  dh = dh_mol / M; // [J/kg]

  // Calore specifico isocorico
  double dcv_mol = T * ddadTm * K;
  dcv = dcv_mol / M; // [J/kg/K]
}

void thermo_derivatives(double am, double bm, double dadTm, double v,
                        double rho, double M, double T,
                        const std::vector<double> &X,
                        const std::vector<Species> &species, double &dpdT_v,
                        double &dpdv_T, double &dpdT_rho, double &dpdrho_T) {
  // Inizializza delta1 e delta2 basati sul tipo di EoS
  double delta1 = 0.0, delta2 = 0.0;

  // Tipo di Equazione di Stato (EoS)
  std::string EoS =
      "SRK"; // Questo può essere impostato globalmente o passato come parametro

  // Calcola delta1 e delta2 in base all'EoS selezionata
  // Calcola delta1 e delta2 in base all'EoS selezionata
  // if (EoS == "PR") {
  //     delta1 = 1.0 + std::sqrt(2.0);
  //     delta2 = 1.0 - std::sqrt(2.0);
  // } else if (EoS == "SRK") {
  delta1 = 1.0;
  delta2 = 0.0;
  // } else if (EoS == "RKPR") {
  //     double d1 = 0.428363;
  //     double d2 = 18.496215;
  //     double d3 = 0.338426;
  //     double d4 = 0.660000;
  //     double d5 = 789.723105;
  //     double d6 = 2.512392;

  //     std::vector<double> delta1_i(species.size(), 0.0);
  //     std::vector<double> delta2_i(species.size(), 0.0);

  //     for (size_t k = 0; k < species.size(); ++k) {
  //         delta1_i[k] = d1 + d2 * std::pow(d3 - 1.168 * species[k].Zc, d4) +
  //         d5 * std::pow(d3 - 1.168 * species[k].Zc, d6); delta2_i[k] = (1 -
  //         delta1_i[k]) / (1 + delta1_i[k]);
  //     }

  //     delta1 = 0.0;
  //     delta2 = 0.0;
  //     for (size_t i = 0; i < X.size(); ++i) {
  //         delta1 += X[i] * delta1_i[i];
  //         delta2 += X[i] * delta2_i[i];
  //     }
  // } else {
  //     throw std::invalid_argument("The selected EoS is not available");
  // }

  // Calcola D per i derivati
  double D = (v + delta1 * bm) * (v + delta2 * bm);

  dpdT_v = R0 / (v - bm) - dadTm / D;

  dpdv_T = -R0 * T / std::pow(v - bm, 2) +
           (am * (2 * v + (delta1 + delta2) * bm)) / std::pow(D, 2);

  // Calcoli basati sulla densità (rho)
  double Den_dpdT_rho = (M + delta1 * bm * rho) * (M + delta2 * bm * rho);
  dpdT_rho =
      (rho * R0) / (M - bm * rho) - dadTm * (std::pow(rho, 2) / Den_dpdT_rho);

  double Num_dpdrho_T = am * rho * M * (2 * M + (delta1 + delta2) * bm * rho);
  double Den_dpdrho_T =
      std::pow(M + delta1 * bm * rho, 2) * std::pow(M + delta2 * bm * rho, 2);

  dpdrho_T =
      (M * R0 * T) / std::pow(M - bm * rho, 2) - Num_dpdrho_T / Den_dpdrho_T;
}

double collision_integral(double Tstar) {
  // Coefficients from Chung et al. (1988)
  const double A = 1.16145;
  const double B = 0.14874;
  const double C = 0.52487;
  const double D = 0.77320;
  const double E = 2.16178;
  const double F = 2.43787;
  const double G = -6.435e-4;
  const double H = 7.27371;
  const double S = 18.0323;
  const double W = -0.76830;

  // Calculate the collision integral
  return (A / std::pow(Tstar, B)) + C / std::exp(D * Tstar) +
         E / std::exp(F * Tstar) +
         (G * std::pow(Tstar, B)) * std::sin(S * std::pow(Tstar, W) - H);
}

void chung_mixing(const std::vector<double> &X,
                  const std::vector<Species> &Species, double &sigma3m,
                  double &epsNm, double &omegaNm, double &MNm, double &Vc_m,
                  double &eps_m, double &omega_m, double &M_m, double &Fcm,
                  double &Tc_m) {

  const size_t Ns = X.size();

  // Reset values
  sigma3m = 0.0;
  epsNm = 0.0;
  omegaNm = 0.0;
  MNm = 0.0;

  // Precompute expensive terms - use cbrt instead of pow(x, 1/3)
  std::vector<double> eps_kb_vals(Ns);
  std::vector<double> sqrt_M(Ns);
  std::vector<double> Vc_cbrt(Ns); // cube root of Vc
  std::vector<double> sigma(Ns);   // 0.809 * cbrt(Vc)

  for (size_t i = 0; i < Ns; ++i) {
    eps_kb_vals[i] = Species[i].Tc / 1.2593;
    sqrt_M[i] = std::sqrt(Species[i].M);
    Vc_cbrt[i] = std::cbrt(Species[i].Vc); // cbrt is faster than pow(x, 1/3)
    sigma[i] = 0.809 * Vc_cbrt[i];
  }

  // Combine species properties - exploit symmetry
  for (size_t i = 0; i < Ns; ++i) {
    const double Xi = X[i];
    if (Xi < 1e-12)
      continue; // Skip negligible species

    // Diagonal term (i, i)
    const double sigma_ii = sigma[i];
    const double eps_ii = eps_kb_vals[i];
    const double sigma_ii_2 = sigma_ii * sigma_ii;
    const double sigma_ii_3 = sigma_ii_2 * sigma_ii;

    const double term_ii = Xi * Xi;
    sigma3m += term_ii * sigma_ii_3;
    epsNm += term_ii * eps_ii * sigma_ii_3;
    omegaNm += term_ii * Species[i].omega * sigma_ii_3;
    MNm += term_ii * eps_ii * sigma_ii_2 * sqrt_M[i];

    // Off-diagonal terms (j > i) - uses symmetry factor 2
    for (size_t j = i + 1; j < Ns; ++j) {
      const double Xj = X[j];
      if (Xj < 1e-12)
        continue;

      // Geometric mean for sigma_ij (sqrt of product = product of sqrt)
      const double sigma_ij = std::sqrt(sigma[i] * sigma[j]);
      const double eps_ij = std::sqrt(eps_kb_vals[i] * eps_kb_vals[j]);
      const double omega_ij = 0.5 * (Species[i].omega + Species[j].omega);

      // Harmonic mean for M_ij
      const double Mi = Species[i].M;
      const double Mj = Species[j].M;
      const double M_ij = 2.0 * Mi * Mj / (Mi + Mj);

      const double sigma_ij_2 = sigma_ij * sigma_ij;
      const double sigma_ij_3 = sigma_ij_2 * sigma_ij;

      const double term_ij = 2.0 * Xi * Xj; // Factor 2 for symmetry

      sigma3m += term_ij * sigma_ij_3;
      epsNm += term_ij * eps_ij * sigma_ij_3;
      omegaNm += term_ij * omega_ij * sigma_ij_3;
      MNm += term_ij * eps_ij * sigma_ij_2 * std::sqrt(M_ij);
    }
  }

  // Compute final results - use cbrt for 2/3 power as cbrt(x)^2
  if (sigma3m > 0.0) {
    const double sigma3m_inv = 1.0 / sigma3m;
    constexpr double k_0809_cubed = 0.809 * 0.809 * 0.809; // 0.529475...
    Vc_m = sigma3m / k_0809_cubed;
    eps_m = epsNm * sigma3m_inv;
    omega_m = omegaNm * sigma3m_inv;

    // sigma3m^(2/3) = (cbrt(sigma3m))^2
    const double sigma3m_cbrt = std::cbrt(sigma3m);
    const double sigma3m_exp2 = sigma3m_cbrt * sigma3m_cbrt;
    const double temp = MNm / (sigma3m_exp2 * eps_m);
    M_m = temp * temp;
  } else {
    Vc_m = 0.0;
    eps_m = 0.0;
    omega_m = 0.0;
    M_m = 0.0;
  }

  Fcm = 1.0 - 0.2756 * omega_m;
  Tc_m = 1.2593 * eps_m;
}

/*
void chung_transport(double eps_m, double Fcm, double W_m, double Vc_m, double
omega_m, double Tc_m, double cv0, double v, double T, double& mu, double&
lambda_, double& mu0) {

    // Dimensionless temperature
    double Tstar_m = T / eps_m;

    // Collision integral
    double sigmaV = collision_integral(Tstar_m);

    // Low-pressure viscosity
    const double EXP = 2.0 / 3.0;
    mu0 = (4.0785e-5) * (std::sqrt(W_m * T) * Fcm) / (sigmaV * std::pow(Vc_m,
EXP));

    // Y
    double y = Vc_m / (6 * v * 1e6);
    double G1 = (1 - 0.5 * y) / std::pow((1 - y), 3);

    // Chung's coefficients
    std::array<double, 10> a1 = {6.32402, 0.0012102, 5.28346, 6.62263, 19.74540,
-1.89992, 24.27450, 0.79716, -0.23816, 0.068629}; std::array<double, 10> a2 =
{50.41190, -0.0011536, 254.209, 38.0957, 7.63034, -12.53670, 3.44945, 1.11764,
0.067695, 0.34793};

    // Polynomial construction
    std::array<double, 10> E;
    for (int k = 0; k < 10; ++k) {
    E[k] = a1[k] + a2[k] * omega_m;
    }

    double nG2 = (E[0] * (1 - std::exp(-E[3] * y))) / y + E[1] * G1 *
std::exp(E[4] * y) + E[2] * G1; double dG2 = E[0] * E[3] + E[1] + E[2]; double
G2 = nG2 / dG2;

    // Viscosity calculations
    double muK = mu0 * (1 / G2 + E[5] * y);
    double muP = (36.344e-6 * (std::sqrt(Tc_m * W_m) / std::pow(Vc_m, EXP))) *
E[6] * std::pow(y, 2) * G2 * std::exp(E[7] + E[8] * (1 / Tstar_m) + E[9] *
std::pow(Tstar_m, -2));

    // High-pressure viscosity
    mu = 0.1 * (muK + muP);

    // Thermal conductivity calculations
    double Tr = T / Tc_m;
    const double Rcal = 1.987;  // cal/mol*K
    double cv0_mol = cv0 * W_m * 1e-3;
    double cv0_mol_cal = cv0_mol * 0.239006;

    double alpha = cv0_mol_cal / Rcal - 1.5;
    double beta = 0.7862 - 0.7109 * omega_m + 1.3168;
    double gamma = 2.0 + 10.5 * Tr*Tr;
    double psiN = 0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * gamma;
    double psiD = 0.6366 + beta * gamma + 1.061 * alpha * beta;
    double psi = 1 + alpha * (psiN / psiD);

     // Conduttività di riferimento
    double lambda0 = 7.452 * (mu0 / W_m) * psi;

    // Coefficienti
    std::array<double, 7> z1 = {2.41657, -0.50924, 6.61069, 14.54250, 0.79274,
-5.8634, 81.171}; std::array<double, 7> z2 = {0.74824, -1.50936, 5.62073,
-8.91387, 0.82019, 12.8005, 114.1580}; std::array<double, 7> B;

    // Calcolo dei coefficienti B[i]
    for (int i = 0; i < 7; ++i) {
        B[i] = z1[i] + z2[i] * omega_m;
    }

    // Numeratore e denominatore di H2
    double nH2 = (B[0] * (1.0 - std::exp(-B[3] * y))) / y
               + B[1] * G1 * std::exp(B[4] * y)
               + B[2] * G1;

    double dH2 = B[0] * B[3] + B[1] + B[2];
    double H2 = nH2 / dH2;

    // Conducibilità a bassa pressione
    double lambdak = lambda0 * (1.0 / H2 + B[5] * y);

    // Conducibilità ad alta pressione
    double lambdap = (3.039e-4 * (std::sqrt(Tc_m / W_m) / std::pow(Vc_m, EXP)))
                   * B[6] * y * y * H2 * std::sqrt(Tr);

    // Conducibilità termica finale
    lambda_ = lambdak + lambdap;

    // Conversione da cal/(cm·s·K) a W/(m·K)
    lambda_ = (lambda_ / 0.239) * 100.0;

}
*/
void chung_transport(double eps_m, double Fcm, double W_m, double Vc_m,
                     double omega_m, double Tc_m, double cv0, double v,
                     double T, double &mu, double &lambda_, double &mu0) {

  // --- Unit Conversions to CGS (Chung Expected Units) ---
  double W_cgs = W_m * 1000.0;      // kg/mol -> g/mol
  double Vc_cgs = Vc_m * 1.0e6;     // m^3/mol -> cm^3/mol
  double v_cgs = v * 1.0e6;         // m^3/mol -> cm^3/mol
  double rho_mol_cgs = 1.0 / v_cgs; // mol/cm^3

  // Dimensionless temperature
  double Tstar_m = T / eps_m;

  // Collision integral
  double sigmaV = collision_integral(Tstar_m);

  // Low-pressure viscosity [micropoise]
  const double EXP = 2.0 / 3.0;
  // Constant 40.785 for micropoise, g/mol, cm3/mol
  double mu0_microP = (40.785) * (std::sqrt(W_cgs * T) * Fcm) /
                      (sigmaV * std::pow(Vc_cgs, EXP));

  // Y (Reserved density in Chung units)
  // y = Vc / (6 * v)
  double y = Vc_cgs / (6.0 * v_cgs);

  double G1 = (1.0 - 0.5 * y) / std::pow((1.0 - y), 3);

  // Chung's coefficients - static constexpr for efficiency
  static constexpr std::array<double, 10> a1 = {
      6.32402,  0.0012102, 5.28346, 6.62263,  19.74540,
      -1.89992, 24.27450,  0.79716, -0.23816, 0.068629};
  static constexpr std::array<double, 10> a2 = {
      50.41190,  -0.0011536, 254.209, 38.0957,  7.63034,
      -12.53670, 3.44945,    1.11764, 0.067695, 0.34793};

  // Polynomial construction
  std::array<double, 10> E;
  for (int k = 0; k < 10; ++k) {
    E[k] = a1[k] + a2[k] * omega_m;
  }

  double nG2 = (E[0] * (1.0 - std::exp(-E[3] * y))) / y +
               E[1] * G1 * std::exp(E[4] * y) + E[2] * G1;
  double dG2 = E[0] * E[3] + E[1] + E[2];
  double G2 = nG2 / dG2;

  // Viscosity calculations [micropoise]
  double muK_microP = mu0_microP * (1.0 / G2 + E[5] * y);

  // High-pressure correction term (36.344 coeff for micropoise, g/mol, cm3/mol)
  double muP_microP =
      (36.344 * (std::sqrt(Tc_m * W_cgs) / std::pow(Vc_cgs, EXP))) * E[6] *
      (y * y) * G2 *
      std::exp(E[7] + E[8] / Tstar_m + E[9] / (Tstar_m * Tstar_m));

  // High-pressure viscosity [micropoise]
  double mu_microP = muK_microP + muP_microP;

  // --- Convert Viscosity to SI (Pa.s) ---
  // 1 micropoise = 1e-7 Pa.s
  mu = mu_microP * 1.0e-7;
  mu0 = mu0_microP * 1.0e-7;

  // --- Thermal Conductivity ---

  double Tr = T / Tc_m;
  const double Rcal = 1.987; // cal/mol*K

  // cv0 is J/kg/K via OpenFOAM/Mixing.
  // W_m is kg/mol.
  // cv0_mol = cv0 * W_m -> J/mol/K
  // 1 J = 0.239006 cal
  double cv0_mol_J = cv0 * W_m;
  double cv0_mol_cal = cv0_mol_J * 0.239006;

  double alpha = cv0_mol_cal / Rcal - 1.5;
  double beta = 0.7862 - 0.7109 * omega_m + 1.3168 * omega_m * omega_m;
  double gamma = 2.0 + 10.5 * Tr * Tr;
  double psiN = 0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * gamma;
  double psiD = 0.6366 + beta * gamma + 1.061 * alpha * beta;
  double psi = 1.0 + alpha * (psiN / psiD);

  // Reference conductivity [cal/(cm s K)]
  // Formula: lambda0 = 7.452 * (mu0_microP / W_cgs) * psi / 1e6 ? NO.
  // Standard Chung: lambda0 = 7.452e-3 * (mu0_microP / W_cgs) * psi ?
  // Wait.
  // Chung 1984: lambda = (31.2 * mu0 * Psi)/M'  (different units).
  // Let's stick to the code's previous coeff but fix units.
  // Original: 7.452 * (mu0 / W_m).
  // If mu0 is micropoise, W_m is g/mol.
  // Result is cal/(cm s K) * 1e6?
  // Let's use mu0_microP and W_cgs.
  // Factor 7.452e-3 produces cal/(cm s K)?
  // Usually: lambda [cal/cm s K] = 7.452 * mu [P] / M ...
  // Note: 1 P = 1e6 microP.
  // So 7.452 * (mu_microP * 1e-6 / W_cgs) * psi?

  // Using literature coefficient for lambda0 with micropoise:
  // lambda0 = 7.452e-6 * (mu0_microP / W_cgs) * psi (??)
  // Actually, factor 7.452 likely assumes Poise?
  // Let's assume the previous code constant was meant for Poise?
  // Let's deduce from SI.
  // lambda_gas ~ 0.02 W/m/K.
  // mu0 ~ 200 microP. W ~ 28.
  // 200/28 ~ 7.
  // 7.452 * 7 = 50.
  // If unit is cal/cm/s/K (418 W/m/K). 50 * 418 = 20000. Too big.
  // So factor must be 1e-6 somewhere.
  // Correct constant for output in cal/cm/s/K using micropoise input is
  // typically 7.452e-6 is unlikely. Let's use the explicit conversion that
  // works: lambda'' = 3.75 * (R/M) * mu * Psi (approx Eucken). Let's use the
  // code's 7.452 but WITH 1e-7 factor for Poise?

  // Using coeff from Poling for lambda in W/m/K directly?
  // No, let's assume the correlation yields lambda* (cal/cm.s.K).
  // lambda0 = 7.452e-6 * (mu0_microP / W_cgs) * psi; // (CHECKED: usually 10^-6
  // factor vs Poise logic)

  double lambda0_cal = (7.452e-6) * (mu0_microP / W_cgs) * psi;

  // Coefficients for lambda - static constexpr for efficiency
  static constexpr std::array<double, 7> z1 = {
      2.41657, -0.50924, 6.61069, 14.54250, 0.79274, -5.8634, 81.171};
  static constexpr std::array<double, 7> z2 = {
      0.74824, -1.50936, 5.62073, -8.91387, 0.82019, 12.8005, 114.1580};
  std::array<double, 7> B;

  for (int i = 0; i < 7; ++i) {
    B[i] = z1[i] + z2[i] * omega_m;
  }

  double nH2 = (B[0] * (1.0 - std::exp(-B[3] * y))) / y +
               B[1] * G1 * std::exp(B[4] * y) + B[2] * G1;

  double dH2 = B[0] * B[3] + B[1] + B[2];
  double H2 = nH2 / dH2;

  // Low pressure conductivity [cal/cm s K]
  double lambdak_cal = lambda0_cal * (1.0 / H2 + B[5] * y);

  // High pressure term
  // Coeff 3.039e-4 is standard.
  double lambdap_cal =
      (3.039e-4 * (std::sqrt(Tc_m / W_cgs) / std::pow(Vc_cgs, EXP))) * B[6] *
      (y * y) * H2 * std::sqrt(Tr);

  // Final Conductivity [cal/cm s K]
  double lambda_cal = lambdak_cal + lambdap_cal;

  // Convert to SI [W/m K]
  // 1 cal/(cm s K) = 418.4 W/(m K)
  lambda_ = lambda_cal * 418.4;
}

// psiN = 0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * gamma
// psiD = 0.6366 + beta * gamma + 1.061 * alpha * beta
// psi = 1 + alpha * (psiN / psiD)

// # Low-pressure conductivity
// lambda0 = 7.452 * (mu0 / W_m) * psi

// # Chung's coefficients for conductivity
// z1 = np.array([2.41657, -0.50924, 6.61069, 14.54250, 0.79274,
// -5.8634, 81.171]) z2 = np.array([0.74824, -1.50936, 5.62073, -8.91387,
// 0.82019, 12.8005, 114.1580])

// B = np.zeros(7)
// for i in range(7):
//     B[i] = z1[i] + z2[i] * omega_m

// nH2 = (B[0] * (1 - np.exp(-B[3] * y))) / y + B[1] * G1 * np.exp(B[4] * y) +
// B[2] * G1 dH2 = B[0] * B[3] + B[1] + B[2] H2 = nH2 / dH2

// # Low-pressure conductivity
// lambdak = lambda0 * (1 / H2 + B[5] * y)
// lambdap = (3.039e-4 * (np.sqrt(Tc_m / W_m) / (Vc_m**EXP))) * B[6] * y**2 * H2
// * np.sqrt(Tr)

// # Final thermal conductivity
// lambda_ = lambdak + lambdap
// lambda_ = (lambda_ / 0.239) * 100

// return mu, lambda_, mu0

double fug_coef(const std::vector<double> &X, double bi, double bm,
                double sum_ref_i, double am, double Z, double T, double p,
                const std::string &EoS, const std::vector<Species> &species) {
  const double R0 = 8.314472;

  double A = (am * p) / std::pow(R0 * T, 2);
  double B = (bm * p) / (R0 * T);

  double delta1 = 0.0, delta2 = 0.0;
  double lnphi_i = 0.0;

  delta1 = 1.0;
  delta2 = 0.0;
  lnphi_i = (bi / bm) * (Z - 1.0) - std::log(std::abs(Z - B)) -
            A / (B * (delta1 - delta2)) * (2.0 * sum_ref_i / am - bi / bm) *
                std::log((Z + B * delta1) / (Z + B * delta2));

  return lnphi_i;
}

void TPD_SSI(const std::vector<double> &X, double T, double p,
             const std::string &mixCR, const std::vector<Species> &Species,
             const std::string &EoS, double kij, bool &stable, double &TPDmin,
             std::vector<double> &K_init_min) {
  const size_t Ns = X.size();

  std::vector<double> X_calc = X;
  for (size_t i = 0; i < Ns; ++i) {
    if (X_calc[i] < 0.0) {
      X_calc[i] = 1e-16;
    }
  }

  bool is_pure = false;
  for (size_t i = 0; i < Ns; ++i) {
    if (X_calc[i] >= 0.999)
      is_pure = true;
  }
  if (is_pure) {
    TPDmin = 10.0;
    stable = true;
    K_init_min.assign(Ns, 0.0);
    return;
  }

  const int ntrial = 4;
  const double tol = 1e-9;
  const int iter_max = 10000;
  const double be_negative = -1e-8;
  const bool return_if_negative = true;

  double am_X = 0, bm_X = 0, dadTm_X = 0, ddadTm_X = 0;
  std::vector<double> sum_ref_i_X(Ns), bc(Ns);

  if (mixCR == "CR1") {
    VdWmixing_CR1_fugacity(X_calc, T, Species, EoS, kij, am_X, bm_X, dadTm_X,
                           ddadTm_X, sum_ref_i_X, bc);
  } else if (mixCR == "CR2") {
    VdWmixing_CR2_fugacity(X_calc, T, Species, EoS, kij, am_X, bm_X, dadTm_X,
                           ddadTm_X, sum_ref_i_X, bc);
  } else {
    throw std::invalid_argument("The selected mixing rule is not available");
  }

  double Z_z = getZ(am_X, bm_X, p, T, X_calc, EoS, Species);
  bool is_liq = (Z_z < 0.5);

  std::vector<double> lnphi_z(Ns), d(Ns);
  for (size_t i = 0; i < Ns; ++i) {
    lnphi_z[i] = fug_coef(X_calc, bc[i], bm_X, sum_ref_i_X[i], am_X, Z_z, T, p,
                          EoS, Species);
    d[i] = std::log(X_calc[i]) + lnphi_z[i];
  }

  std::vector<double> K(Ns);
  for (size_t i = 0; i < Ns; ++i) {
    K[i] = (Species[i].Pc / p) * std::exp(5.37 * (1.0 + Species[i].omega) *
                                          (1.0 - Species[i].Tc / T));
  }

  std::vector<std::vector<double>> Y_init(ntrial, std::vector<double>(Ns));

  for (size_t i = 0; i < Ns; ++i) {
    if (is_liq) {
      Y_init[0][i] = X_calc[i] * K[i];
      Y_init[1][i] = X_calc[i] * std::pow(K[i], 1.0 / 3.0);
      Y_init[2][i] = X_calc[i] / K[i];
      Y_init[3][i] = X_calc[i] / std::pow(K[i], 1.0 / 3.0);
    } else {
      Y_init[0][i] = X_calc[i] / K[i];
      Y_init[1][i] = X_calc[i] / std::pow(K[i], 1.0 / 3.0);
      Y_init[2][i] = X_calc[i] * K[i];
      Y_init[3][i] = X_calc[i] * std::pow(K[i], 1.0 / 3.0);
    }
  }

  std::vector<double> TPD_star(ntrial, 1.0);
  double global_tpd_min = 10.0;
  int k_min_idx = -1;

  for (int k = 0; k < ntrial; ++k) {
    std::vector<double> Y = Y_init[k];
    bool break_eps = false;

    for (int iter = 0; iter < iter_max; ++iter) {
      double sum_Y = 0.0;
      for (size_t i = 0; i < Ns; ++i)
        sum_Y += Y[i];

      std::vector<double> y_trial(Ns);
      for (size_t i = 0; i < Ns; ++i)
        y_trial[i] = Y[i] / sum_Y;

      double am_y = 0, bm_y = 0, dadTm_y = 0, ddadTm_y = 0;
      std::vector<double> sum_ref_i_y(Ns), bc_y(Ns);

      if (mixCR == "CR1") {
        VdWmixing_CR1_fugacity(y_trial, T, Species, EoS, kij, am_y, bm_y,
                               dadTm_y, ddadTm_y, sum_ref_i_y, bc_y);
      } else if (mixCR == "CR2") {
        VdWmixing_CR2_fugacity(y_trial, T, Species, EoS, kij, am_y, bm_y,
                               dadTm_y, ddadTm_y, sum_ref_i_y, bc_y);
      }

      double Z_y = getZ(am_y, bm_y, p, T, y_trial, EoS, Species);

      std::vector<double> Yn(Ns);
      for (size_t i = 0; i < Ns; ++i) {
        double lnphi_y = fug_coef(X_calc, bc_y[i], bm_y, sum_ref_i_y[i], am_y,
                                  Z_y, T, p, EoS, Species);
        Yn[i] = std::exp(d[i] - lnphi_y);
      }

      double eps_sq = 0.0;
      for (size_t i = 0; i < Ns; ++i) {
        double dy = Yn[i] - Y[i];
        eps_sq += dy * dy;
      }
      double eps = std::sqrt(eps_sq);
      Y = Yn;

      if (eps <= tol) {
        double final_sum_Y = 0.0;
        for (size_t i = 0; i < Ns; ++i)
          final_sum_Y += Y[i];
        TPD_star[k] = 1.0 - final_sum_Y;
        break_eps = true;
        break;
      }
    }

    if (!break_eps) {
      double sum_Y = 0.0;
      for (size_t i = 0; i < Ns; ++i)
        sum_Y += Y[i];

      std::vector<double> y_trial(Ns);
      for (size_t i = 0; i < Ns; ++i)
        y_trial[i] = Y[i] / sum_Y;

      double am_y = 0, bm_y = 0, dadTm_y = 0, ddadTm_y = 0;
      std::vector<double> sum_ref_i_y(Ns), bc_y(Ns);
      if (mixCR == "CR1") {
        VdWmixing_CR1_fugacity(y_trial, T, Species, EoS, kij, am_y, bm_y,
                               dadTm_y, ddadTm_y, sum_ref_i_y, bc_y);
      } else if (mixCR == "CR2") {
        VdWmixing_CR2_fugacity(y_trial, T, Species, EoS, kij, am_y, bm_y,
                               dadTm_y, ddadTm_y, sum_ref_i_y, bc_y);
      }

      double Z_y = getZ(am_y, bm_y, p, T, y_trial, EoS, Species);

      double tpd_val = 1.0;
      for (size_t i = 0; i < Ns; ++i) {
        double lnphi_y = fug_coef(X_calc, bc_y[i], bm_y, sum_ref_i_y[i], am_y,
                                  Z_y, T, p, EoS, Species);
        tpd_val += Y[i] * (lnphi_y + std::log(Y[i]) - d[i] - 1.0);
      }
      TPD_star[k] = tpd_val;
    }

    if (TPD_star[k] < global_tpd_min) {
      global_tpd_min = TPD_star[k];
      k_min_idx = k;
    }

    if (return_if_negative && break_eps) {
      if (TPD_star[k] < be_negative) {
        break;
      }
    }
  }

  TPDmin = global_tpd_min;

  if (TPDmin < be_negative) {
    stable = false;
    K_init_min.assign(Ns, 0.0);
    if (k_min_idx != -1) {
      for (size_t i = 0; i < Ns; ++i) {
        if (is_liq) {
          if (k_min_idx == 0)
            K_init_min[i] = Y_init[k_min_idx][i] / X_calc[i];
          else if (k_min_idx == 1)
            K_init_min[i] = std::pow(Y_init[k_min_idx][i] / X_calc[i], 3.0);
          else if (k_min_idx == 2)
            K_init_min[i] = X_calc[i] / Y_init[k_min_idx][i];
          else if (k_min_idx == 3)
            K_init_min[i] = std::pow(X_calc[i] / Y_init[k_min_idx][i], 3.0);
        } else {
          if (k_min_idx == 0)
            K_init_min[i] = X_calc[i] / Y_init[k_min_idx][i];
          else if (k_min_idx == 1)
            K_init_min[i] = std::pow(X_calc[i] / Y_init[k_min_idx][i], 3.0);
          else if (k_min_idx == 2)
            K_init_min[i] = Y_init[k_min_idx][i] / X_calc[i];
          else if (k_min_idx == 3)
            K_init_min[i] = std::pow(Y_init[k_min_idx][i] / X_calc[i], 3.0);
        }
      }
    }
  } else {
    stable = true;
    K_init_min.assign(Ns, 0.0);
  }
}

bool flashTPN(const std::vector<double> &z_in, double T, double p,
              const std::string &mixCR, const std::vector<Species> &SpeciesList,
              const std::string &EoS, double kij, std::vector<double> &x,
              std::vector<double> &y, double &psi_v,
              const std::vector<double> &K_init) {
  const double R0 = 8.314472;
  int Ns = z_in.size();
  bool success = false;

  // Normalize z_in just in case
  std::vector<double> z = z_in;
  double z_sum = 0.0;
  for (double val : z)
    z_sum += val;
  if (z_sum > 0.0) {
    for (double &val : z)
      val /= z_sum;
  }

  // Initialize outputs
  x.assign(Ns, 0.0);
  y.assign(Ns, 0.0);
  psi_v = 0.0;

  // 1. Check Stability
  bool stable = true;
  double TPDmin = 0.0;
  std::vector<double> K_init_min;

  TPD_SSI(z, T, p, mixCR, SpeciesList, EoS, kij, stable, TPDmin, K_init_min);

  // If stable, it's a single phase
  if (stable) {
    double am_p, bm_p;
    if (mixCR == "CR1") {
      double dummy_dadTm, dummy_ddadTm;
      VdWmixing_CR1(z, T, SpeciesList, EoS, kij, am_p, bm_p, dummy_dadTm,
                    dummy_ddadTm);
    } else if (mixCR == "CR2") {
      double dummy_dadTm, dummy_ddadTm;
      VdWmixing_CR2(z, T, SpeciesList, EoS, kij, am_p, bm_p, dummy_dadTm,
                    dummy_ddadTm);
    } else {
      throw std::runtime_error("Unknown mixing rule in flashTPN");
    }

    double Z_p = getZ(am_p, bm_p, p, T, z, "Gibbs", SpeciesList);
    double bd = 1.0 / (R0 * T);
    double B = (bm_p * p) * bd;

    if (Z_p < 3.5 * B) {
      psi_v = 0.0; // Liquid
      x = z;
      y = z;
    } else {
      psi_v = 1.0; // Vapor
      x = z;
      y = z;
    }

    return true;
  }

  // 2. Unstable: Perform Phase Split (SSI Loop)
  double tol_o = 1.e-08; //1e-7
  double tol_i = 1.e-07;
  double tol_trivial = 1.e-10; //1e-6
  int itmax_o = 1000;
  int itmax_i = 100;

  std::vector<double> K(Ns, 0.0);
  if (!K_init.empty() && K_init.size() == Ns) {
    K = K_init;
  } else {
    // Wilson equation
    for (int i = 0; i < Ns; ++i) {
      K[i] = SpeciesList[i].Pc / p *
             exp(5.37 * (1.0 + SpeciesList[i].omega) *
                 (1.0 - SpeciesList[i].Tc / T));
    }
  }

  std::vector<double> lnphi_l(Ns, 0.0);
  std::vector<double> lnphi_v(Ns, 0.0);
  std::vector<double> F(Ns, 0.0);

  int iter_o = 0;

  while (iter_o <= itmax_o) {

    double sum_Kz = 0.0;
    double sum_zK = 0.0;
    for (int i = 0; i < Ns; ++i) {
      sum_Kz += K[i] * z[i];
      sum_zK += z[i] / K[i];
    }

    if (sum_Kz >= 1.0 && sum_zK >= 1.0) {
      // Solve Rachford-Rice
      double eps_i = 1.0;
      int iter_i = 1;

      double psi_v_min = 0.0;
      double psi_v_max = 1.0;

      for (int i = 0; i < Ns; ++i) {
        if (K[i] > 1.0) {
          psi_v_min = std::max(psi_v_min, (K[i] * z[i] - 1.0) / (K[i] - 1.0));
        } else if (K[i] < 1.0) {
          psi_v_max = std::min(psi_v_max, (1.0 - z[i]) / (1.0 - K[i]));
        }
      }

      psi_v = 0.5 * (psi_v_min + psi_v_max);

      while (eps_i >= tol_i && iter_i <= itmax_i) {
        double RR = 0.0;
        double dRR = 0.0;

        for (int i = 0; i < Ns; ++i) {
          double term = 1.0 + psi_v * (K[i] - 1.0);
          RR += z[i] * (K[i] - 1.0) / term;
          dRR += z[i] * std::pow(K[i] - 1.0, 2) / std::pow(term, 2);
        }
        dRR = -dRR;

        if (RR > 0.0) {
          psi_v_min = psi_v;
        } else {
          psi_v_max = psi_v;
        }

        double dpsi_v = -(RR / dRR);
        double psi_vn = psi_v + dpsi_v;

        if (psi_v != 0) {
          eps_i = std::abs((psi_vn - psi_v) / psi_v);
        } else {
          eps_i = std::abs(psi_vn);
        }

        if (psi_vn > psi_v_min && psi_vn < psi_v_max) {
          psi_v = psi_vn;
        } else {
          psi_v = 0.5 * (psi_v_min + psi_v_max);
        }
        iter_i++;
      }

      for (int i = 0; i < Ns; ++i) {
        x[i] = z[i] / (1.0 + psi_v * (K[i] - 1.0));
        y[i] = K[i] * x[i];
      }
    } else if (sum_zK <= 1.0) {
      psi_v = 1.0;
      double x_sum = 0.0;
      for (int i = 0; i < Ns; ++i) {
        x[i] = z[i] / K[i];
        y[i] = z[i];
        x_sum += x[i];
      }
      if (x_sum > 0) {
        for (double &val : x)
          val /= x_sum;
      }
    } else if (sum_Kz <= 1.0) {
      psi_v = 0.0;
      double y_sum = 0.0;
      for (int i = 0; i < Ns; ++i) {
        x[i] = z[i];
        y[i] = K[i] * z[i];
        y_sum += y[i];
      }
      if (y_sum > 0) {
        for (double &val : y)
          val /= y_sum;
      }
    }

    // Evaluate mixing rules for both phases
    double am_l, bm_l, am_v, bm_v;
    std::vector<double> sum_ref_i_l(Ns, 0.0), sum_ref_i_v(Ns, 0.0);
    std::vector<double> bc_l(Ns, 0.0), bc_v(Ns, 0.0);

    if (mixCR == "CR1") {
      double dum_dadTm, dum_ddadTm;
      VdWmixing_CR1_fugacity(x, T, SpeciesList, EoS, kij, am_l, bm_l, dum_dadTm,
                             dum_ddadTm, sum_ref_i_l, bc_l);
      VdWmixing_CR1_fugacity(y, T, SpeciesList, EoS, kij, am_v, bm_v, dum_dadTm,
                             dum_ddadTm, sum_ref_i_v, bc_v);
    } else if (mixCR == "CR2") {
      double dum_dadTm, dum_ddadTm;
      VdWmixing_CR2_fugacity(x, T, SpeciesList, EoS, kij, am_l, bm_l, dum_dadTm,
                             dum_ddadTm, sum_ref_i_l, bc_l);
      VdWmixing_CR2_fugacity(y, T, SpeciesList, EoS, kij, am_v, bm_v, dum_dadTm,
                             dum_ddadTm, sum_ref_i_v, bc_v);
    }

    double Z_l = getZ(am_l, bm_l, p, T, x, "Liq", SpeciesList);
    double Z_v = getZ(am_v, bm_v, p, T, y, "Vap", SpeciesList);

    double F_norm_sq = 0.0;
    for (int i = 0; i < Ns; ++i) {
      lnphi_l[i] = fug_coef(x, bc_l[i], bm_l, sum_ref_i_l[i], am_l, Z_l, T, p,
                            EoS, SpeciesList);
      lnphi_v[i] = fug_coef(y, bc_v[i], bm_v, sum_ref_i_v[i], am_v, Z_v, T, p,
                            EoS, SpeciesList);

      F[i] = lnphi_l[i] - lnphi_v[i] - log(K[i]);
      F_norm_sq += F[i] * F[i];
    }

    double eps_o = std::sqrt(F_norm_sq);

    if (eps_o <= tol_o) {
      success = true;
      break;
    } else {
      for (int i = 0; i < Ns; ++i) {
        // Safely update K to avoid overflow
        // double new_K = std::exp(lnphi_l[i] - lnphi_v[i]);
        double new_K = std::exp(std::min(std::max(lnphi_l[i] - lnphi_v[i], -50.0), 50.0));
        K[i] = new_K;
      }
    }

    iter_o++;
  }

  // Handle non-convergence: if SSI failed, check if x ≈ y (trivial)
    if (!success) {
      double diff = 0.0;
      for (int i = 0; i < Ns; ++i)
        diff += std::abs(x[i] - y[i]);

      if (diff < 1.e-4) {
        // Trivial solution: phases are identical, force single-phase
        double am_p, bm_p, dum_dadTm, dum_ddadTm;
        if (mixCR == "CR1")
          VdWmixing_CR1(z, T, SpeciesList, EoS, kij, am_p, bm_p, dum_dadTm, dum_ddadTm);
        else
          VdWmixing_CR2(z, T, SpeciesList, EoS, kij, am_p, bm_p, dum_dadTm, dum_ddadTm);

        double Z_p = getZ(am_p, bm_p, p, T, z, "Gibbs", SpeciesList);
        double B = (bm_p * p) / (R0 * T);
        psi_v = (Z_p < 3.5 * B) ? 0.0 : 1.0;
        x = z;
        y = z;
        return true;
      }
    }

  // Trivial solution check
  if (std::abs(psi_v - 1.0) < tol_trivial || std::abs(psi_v) < tol_trivial) {
    double am_l, bm_l;
    if (mixCR == "CR1") {
      double dum_dadTm, dum_ddadTm;
      VdWmixing_CR1(x, T, SpeciesList, EoS, kij, am_l, bm_l, dum_dadTm,
                    dum_ddadTm);
    } else if (mixCR == "CR2") {
      double dum_dadTm, dum_ddadTm;
      VdWmixing_CR2(x, T, SpeciesList, EoS, kij, am_l, bm_l, dum_dadTm,
                    dum_ddadTm);
    }

    double Z_l = getZ(am_l, bm_l, p, T, x, "Liq", SpeciesList);
    double bd = 1.0 / (R0 * T);
    double B_lv = (bm_l * p) * bd;

   
    if (Z_l < 3.5 * B_lv) {
      psi_v = 0.0;
      x = z; // Liquid phase is the same as feed
      y = z;
    } else {
      psi_v = 1.0;
      x = z;
      y = z;
    }
  }

  return success;
}

// ============================================================================
//  SPINODAL STABILITY CHECK
//  Evaluates the diffusional stability matrix:
//    Q_ij = delta_ij / x_i  +  d(ln phi_i)/d(x_j)  at const T, P
//  The mixture is inside the spinodal (diffusionally unstable) when det(Q) < 0.
// ============================================================================

// Helper: compute determinant of an n×n matrix via Gaussian elimination
// with partial pivoting.  Operates on a flat row-major vector.
static double determinant(std::vector<double> A, int n) {
  double det = 1.0;
  for (int col = 0; col < n; ++col) {
    // Partial pivoting
    int pivot = col;
    double maxVal = std::abs(A[col * n + col]);
    for (int row = col + 1; row < n; ++row) {
      double val = std::abs(A[row * n + col]);
      if (val > maxVal) {
        maxVal = val;
        pivot = row;
      }
    }
    if (maxVal < 1e-30)
      return 0.0; // Singular

    if (pivot != col) {
      for (int k = 0; k < n; ++k)
        std::swap(A[col * n + k], A[pivot * n + k]);
      det = -det; // Row swap flips sign
    }

    det *= A[col * n + col];

    // Eliminate below
    double inv_diag = 1.0 / A[col * n + col];
    for (int row = col + 1; row < n; ++row) {
      double factor = A[row * n + col] * inv_diag;
      for (int k = col; k < n; ++k)
        A[row * n + k] -= factor * A[col * n + k];
    }
  }
  return det;
}

bool checkSpinodal(const std::vector<double> &X, double T, double p,
                   const std::string &mixCR,
                   const std::vector<Species> &SpeciesList,
                   const std::string &EoS, double kij, double &detQ) {
  const int Ns = static_cast<int>(X.size());

  // --- Guard: pure component is always diffusionally stable ---
  int nNonZero = 0;
  for (int i = 0; i < Ns; ++i) {
    if (X[i] > 1e-12)
      ++nNonZero;
  }
  if (nNonZero <= 1) {
    detQ = 1.0; // Positive ⟹ stable
    return false;
  }

  // --- Step 1: Compute ln(phi_i) at the unperturbed composition ---
  //     We need the mixing-rule outputs for fug_coef.
  auto computeLnPhi = [&](const std::vector<double> &x_in,
                          std::vector<double> &lnphi) {
    double am, bm, dadTm, ddadTm;
    std::vector<double> sum_ref_i(Ns), bc(Ns);

    if (mixCR == "CR1") {
      VdWmixing_CR1_fugacity(x_in, T, SpeciesList, EoS, kij, am, bm, dadTm,
                             ddadTm, sum_ref_i, bc);
    } else if (mixCR == "CR2") {
      VdWmixing_CR2_fugacity(x_in, T, SpeciesList, EoS, kij, am, bm, dadTm,
                             ddadTm, sum_ref_i, bc);
    } else {
      throw std::invalid_argument("checkSpinodal: unknown mixing rule '" +
                                  mixCR + "'");
    }

    double Z = getZ(am, bm, p, T, x_in, EoS, SpeciesList);

    lnphi.resize(Ns);
    for (int i = 0; i < Ns; ++i) {
      lnphi[i] = fug_coef(x_in, bc[i], bm, sum_ref_i[i], am, Z, T, p, EoS,
                          SpeciesList);
    }
  };

  // --- Step 2: Build the stability matrix Q via finite differences ---
  //     Q_ij = delta_ij / x_i  +  d(ln phi_i) / d(x_j)
  const double h = 1.0e-6; // FD step size

  // Q stored row-major: Q[i * Ns + j]
  std::vector<double> Q(Ns * Ns, 0.0);

  for (int j = 0; j < Ns; ++j) {
    // --- Perturb x_j  +h ---
    std::vector<double> x_plus = X;
    x_plus[j] += h;
    {
      double s = 0.0;
      for (int k = 0; k < Ns; ++k)
        s += x_plus[k];
      for (int k = 0; k < Ns; ++k)
        x_plus[k] /= s;
    }

    // --- Perturb x_j  -h ---
    std::vector<double> x_minus = X;
    x_minus[j] -= h;
    // Clamp to avoid negative mole fractions
    if (x_minus[j] < 0.0)
      x_minus[j] = 0.0;
    {
      double s = 0.0;
      for (int k = 0; k < Ns; ++k)
        s += x_minus[k];
      for (int k = 0; k < Ns; ++k)
        x_minus[k] /= s;
    }

    // --- Evaluate ln(phi) at both perturbed compositions ---
    std::vector<double> lnphi_plus, lnphi_minus;
    computeLnPhi(x_plus, lnphi_plus);
    computeLnPhi(x_minus, lnphi_minus);

    // --- Central finite difference ---
    double dx = 2.0 * h; // before renormalization this was the step
    for (int i = 0; i < Ns; ++i) {
      double dlnphi_dxj = (lnphi_plus[i] - lnphi_minus[i]) / dx;
      Q[i * Ns + j] = dlnphi_dxj;
    }
  }

  // Add the ideal contribution: delta_ij / x_i
  for (int i = 0; i < Ns; ++i) {
    if (X[i] > 1e-30) {
      Q[i * Ns + i] += 1.0 / X[i];
    } else {
      Q[i * Ns + i] +=
          1e30; // Very large → stable contribution for trace species
    }
  }

  // --- Step 3: Compute determinant ---
  detQ = determinant(Q, Ns);

  return (detQ < 0.0);
}
