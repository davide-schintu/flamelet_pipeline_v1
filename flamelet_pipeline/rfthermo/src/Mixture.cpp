#include "../include/Mixture.h"
#include <chrono> // assicurati di includerlo
#include <cmath>
#include <iomanip> // per std::setprecision e std::fixed
#include <limits>

// Costruttore che inizializza la classe Mixture e legge il database delle
// specie
Mixture::Mixture(std::string database_name,
                 const std::vector<std::string> &orderedSpeciesNames) {
  // Legge i dati delle specie dal database
  this->speciesList = database_reader(database_name, orderedSpeciesNames);

  // Crea una lista dei nomi delle specie
  for (const auto &specie : this->speciesList) {
    Species_names.push_back(specie.Names);
  }

  // Inizializza la composizione molare (Y) a zero
  Y.resize(speciesList.size(), 0.0);
  X.resize(speciesList.size(), 0.0);

  // Imposta la temperatura e la pressione a zero
  T = 0.0; // Così se ci sono errori si sgama subito
  P = 0.0;
  R = 0.0;
  EoS = "SRK";
}

// //Aggiungere direttamente output di h, mu, rho ecc.
// double Mixture::solveTemperatureFromH(double h_target,
//                              double T_low = 80.0, double T_high = 5000.0,
//                              double tol = 1e-4, int max_iter = 100) {
//     double T_mid, h_mid;
//     double Z, rho, cp, mu, lambda, alpha;

//     for (int i = 0; i < max_iter; ++i) {
//         T_mid = 0.5 * (T_low + T_high);
//         setTP(T_mid, this->P);
//         calculateProperties(Z, rho, h_mid, cp, mu, lambda, alpha);

//         double err = (h_mid - h_target);

//         if (std::abs(err) < tol) {
//             return T_mid; // Converged
//         }

//         // Decide direction
//         setTP(T_low, this->P);
//         double h_low;
//         calculateProperties(Z, rho, h_low, cp, mu, lambda, alpha);

//         if ((h_low - h_target) * err < 0) {
//             T_high = T_mid;
//         } else {
//             T_low = T_mid;
//         }
//         // std::cout << i << std::endl;
//     }

//     throw std::runtime_error("Temperature solver did not converge.");
// }

double Mixture::solveTemperatureFromH(double h_target, double T_guess,
                                      double tol, int max_iter) {
  static bool logged_marker = false;
  if (!logged_marker) {
    std::cerr << "[RFThermo debug marker 2026-05-11] "
                 "solveTemperatureFromH mono-phase/newton path"
              << std::endl;
    logged_marker = true;
  }

  double T = T_guess;
  if (T < 50.0)
    T = 300.0; // Safety reset for bad guesses

  double Z, rho, h, cp;
  double am = 0.0, bm = 0.0, dadTm = 0.0, ddadTm = 0.0;

  // Pre-calculate bm (constant) if possible, but VdWMixing does all.
  // Optimization: VdWMixing is fast enough for now.

  for (int i = 0; i < max_iter; ++i) {
    setTP(T, this->P);

    // Update derivatives and mixing rules
    VdWMixing(am, bm, dadTm, ddadTm);

    // Calculate Z (Cubic EOS)
    Z = calculateZ(am, bm);

    // Calculate Density
    double v;
    getRho(Z, this->T, this->P, this->Mavg, rho, v);

    // Calculate Ideal Gas Properties
    double cp0, cv0, h0;
    getIdealDepartureProperties(this->T, this->Y, this->speciesList, this->Mavg,
                                cp0, cv0, h0);

    // Calculate Departure Functions
    double de, dh, dcv;
    departureFunctions(am, bm, dadTm, ddadTm, this->Mavg, v, this->T, this->P,
                       this->Y, this->EoS, this->speciesList, de, dh, dcv);

    h = h0 + dh;

    // Convergence Check
    double err = (h - h_target);
    if (std::abs(err) < tol || std::abs(err / (h_target + 1e-10)) < 1e-6) {
      return T; // Converged
    }

    // Calculate Cp for Newton Step
    double cv = cv0 + dcv;
    double dpdT_v, dpdv_T, dpdT_rho, dpdrho_T;
    thermo_derivatives(am, bm, dadTm, v, rho, this->Mavg, this->T, this->X,
                       this->speciesList, dpdT_v, dpdv_T, dpdT_rho, dpdrho_T);
    cp = cv + (T / (rho * rho)) * ((dpdT_rho * dpdT_rho) / dpdrho_T);

    // Newton Step
    if (cp < 1e-6)
      cp = 1.0; // Avoid division by zero
    double dT = -err / cp;

    // Step Limiting (Safety)
    double max_step = 0.2 * T; // Limit change to 20% of T
    if (std::abs(dT) > max_step)
      dT = (dT > 0 ? max_step : -max_step);

    T += dT;

    // Bounds Checck
    if (T < 50.0)
      T = 50.0;
    if (T > 5000.0)
      T = 5000.0;
  }

  // Fallback or Error
  // std::cerr << "Warning: solveTemperatureFromH did not converge fully. Error:
  // " << (h - h_target) << std::endl;
  return T; // Return best estimate
}

double Mixture::solveTemperatureFromH_VLE(double h_target, double T_guess,
                                          double tol, int max_iter, bool& is_converged) {
  static bool logged_marker = false;
  if (!logged_marker) {
    std::cerr << "[RFThermo debug marker 2026-05-11] "
                 "solveTemperatureFromH_VLE VLE-root Liq/Vap path v3"
              << std::endl;
    logged_marker = true;
  }

  double h = 0.0;
  double T = T_guess;
  // if (T < 50.0)
  //   T = 300.0; // Safety reset

  bool converged = false;
  double eps = 1.0;
  int iter = 0;

  double p = this->P;
  std::string mixCR = "CR1"; // make this an argument later if needed
  double kij = 0.0;

  while (iter < max_iter) {
    iter++;

    // 1) Flash problem, (X,T,p)-> (x, y, T, p) and Vapour phase fraction
    std::vector<double> x, y;
    double beta_mol = 0.0;
    std::vector<double> dummy_K_init;
    bool success = flashTPN(this->X, T, p, mixCR, this->speciesList, this->EoS,
                            kij, x, y, beta_mol, dummy_K_init);

    // 2) Mass-based Phase Fraction
    double M_l = 0.0;
    double M_v = 0.0;
    for (size_t i = 0; i < this->X.size(); ++i) {
      M_l += x[i] * this->speciesList[i].M;
      M_v += y[i] * this->speciesList[i].M;
    }

    double beta_mass = 0.0;
    if (this->Mavg > 0) {
      beta_mass = beta_mol * (M_v / this->Mavg);
    }

    // 3) VdW Mixing-Rules
    double am_l, bm_l, dadTm_l, ddadTm_l;
    double am_v, bm_v, dadTm_v, ddadTm_v;

    if (mixCR == "CR1") {
      VdWmixing_CR1(x, T, this->speciesList, this->EoS, kij, am_l, bm_l,
                    dadTm_l, ddadTm_l);
      VdWmixing_CR1(y, T, this->speciesList, this->EoS, kij, am_v, bm_v,
                    dadTm_v, ddadTm_v);
    } else if (mixCR == "CR2") {
      VdWmixing_CR2(x, T, this->speciesList, this->EoS, kij, am_l, bm_l,
                    dadTm_l, ddadTm_l);
      VdWmixing_CR2(y, T, this->speciesList, this->EoS, kij, am_v, bm_v,
                    dadTm_v, ddadTm_v);
    } else {
      throw std::invalid_argument("The selected mixing rule is not available");
    }

    // 4) Compute liquid and vapor compressibility factor
    double Z_l = getZ(am_l, bm_l, p, T, x, "Liq", this->speciesList);
    double Z_v = getZ(am_v, bm_v, p, T, y, "Vap", this->speciesList);

    double v_l = (Z_l * R0 * T) / p;
    double v_v = (Z_v * R0 * T) / p;

    double rho_l = M_l / v_l;
    double rho_v = M_v / v_v;
    // 5) Compute liquid and vapor ideal-gas enthalpy
    std::vector<double> x_mass(this->X.size(), 0.0);
    std::vector<double> y_mass(this->X.size(), 0.0);

    for (size_t i = 0; i < this->X.size(); ++i) {
      if (M_l > 0)
        x_mass[i] = (x[i] * this->speciesList[i].M) / M_l;
      if (M_v > 0)
        y_mass[i] = (y[i] * this->speciesList[i].M) / M_v;
    }

    double cp0_l = 0.0, cv0_l = 0.0, h0_l = 0.0;
    double cp0_v = 0.0, cv0_v = 0.0, h0_v = 0.0;

    if (M_l > 0)
      getIdealDepartureProperties(T, x_mass, this->speciesList, M_l, cp0_l,
                                  cv0_l, h0_l);
    if (M_v > 0)
      getIdealDepartureProperties(T, y_mass, this->speciesList, M_v, cp0_v,
                                  cv0_v, h0_v);

    // cv0 è già calcolato correttamente da getIdealDepartureProperties
    // come cp0 - R0/M (con R0 e M in unità coerenti). NON sovrascrivere.

    // 6) Compute Departure Functions
    double de_l, dh_l = 0.0, dcv_l = 0.0;
    double de_v, dh_v = 0.0, dcv_v = 0.0;

    if (M_l > 0 && v_l > 0) {
      departureFunctions(am_l, bm_l, dadTm_l, ddadTm_l, M_l, v_l, T, p, x_mass,
                         this->EoS, this->speciesList, de_l, dh_l, dcv_l);
    }
    if (M_v > 0 && v_v > 0) {
      departureFunctions(am_v, bm_v, dadTm_v, ddadTm_v, M_v, v_v, T, p, y_mass,
                         this->EoS, this->speciesList, de_v, dh_v, dcv_v);
    }

    // 7) Compute enthalpy in liquid and vapor phase
    double h_l = h0_l + dh_l;
    double h_v = h0_v + dh_v;
    double cv_l = cv0_l + dcv_l;
    double cv_v = cv0_v + dcv_v;

    // Mixture enthalpy (mass basis)
    h = beta_mass * h_v + (1.0 - beta_mass) * h_l;

    // 8) Compute the residuum
    eps = std::abs((h_target - h) / h_target);

    if (eps <= tol) {
      converged = true;
      break;
    }

    // 9) Thermo-Derivatives for Cp
    double dpdT_v_l = 0.0, dpdv_T_l = -1.0, dpdT_rho_l = 0.0, dpdrho_T_l = 1.0;
    double cp_l = 0.0;

    if (M_l > 0 && v_l > 0) {
      thermo_derivatives(am_l, bm_l, dadTm_l, v_l, rho_l, M_l, T, x,
                         this->speciesList, dpdT_v_l, dpdv_T_l, dpdT_rho_l,
                         dpdrho_T_l);
      if (dpdv_T_l != 0) {
        cp_l = cv_l - (T * dpdT_v_l * dpdT_v_l / dpdv_T_l) / M_l;
      } else {
        cp_l = cv_l;
      }
    }

    double dpdT_v_v = 0.0, dpdv_T_v = -1.0, dpdT_rho_v = 0.0, dpdrho_T_v = 1.0;
    double cp_v = 0.0;

    if (M_v > 0 && v_v > 0) {
      thermo_derivatives(am_v, bm_v, dadTm_v, v_v, rho_v, M_v, T, y,
                         this->speciesList, dpdT_v_v, dpdv_T_v, dpdT_rho_v,
                         dpdrho_T_v);
      if (dpdv_T_v != 0) {
        cp_v = cv_v - (T * dpdT_v_v * dpdT_v_v / dpdv_T_v) / M_v;
      } else {
        cp_v = cv_v;
      }
    }

    double cp = beta_mass * cp_v + (1.0 - beta_mass) * cp_l;

    // UPDATE TEMPERATURE
    double alpha = 1.0e-01;
    double delta_h = (h - h_target);
    double dT = -delta_h / cp;
    double C = 1.0 / (1.0 + std::abs(dT * alpha));
    T += dT * C;

    // Bounds Check
    if (T < 50.0)
      T = 50.0;
    if (T > 5000.0)
      T = 5000.0;
    ;
  }

  if (!converged) {
    struct EvalResult {
      bool valid;
      double T;
      double h;
      double residual;
    };

    auto eval_at_T = [&](double T_eval) {
      EvalResult result{false, T_eval, 0.0,
                        std::numeric_limits<double>::quiet_NaN()};
      if (T_eval < 50.0 || T_eval > 5000.0 || !std::isfinite(T_eval)) {
        return result;
      }

      double Z_eval = 0.0, rho_eval = 0.0, h_eval = 0.0, cp_eval = 0.0;
      double mu_eval = 0.0, lambda_eval = 0.0, alpha_eval = 0.0;
      double dpdrho_T_eval = 0.0, sound_speed_eval = 0.0, beta_mass_eval = 0.0;

      setTP(T_eval, p);
      calculateProperties_VLE(Z_eval, rho_eval, h_eval, cp_eval, mu_eval,
                              lambda_eval, alpha_eval, dpdrho_T_eval,
                              sound_speed_eval, beta_mass_eval);

      if (std::isfinite(h_eval)) {
        result.valid = true;
        result.h = h_eval;
        result.residual = h_eval - h_target;
      }
      return result;
    };

    auto rel_error = [&](double residual) {
      return std::abs(residual / (h_target + 1e-10));
    };

    EvalResult best = eval_at_T(T_guess);
    if (best.valid && rel_error(best.residual) <= tol) {
      T = best.T;
      h = best.h;
      converged = true;
    } else {
      double step = std::max(0.5, 0.02 * std::max(T_guess, 50.0));
      EvalResult prev_low = best;
      EvalResult prev_high = best;
      EvalResult bracket_low{false, 0.0, 0.0, 0.0};
      EvalResult bracket_high{false, 0.0, 0.0, 0.0};
      bool bracket_found = false;

      auto update_best = [&](const EvalResult &candidate) {
        if (!candidate.valid) {
          return;
        }
        if (!best.valid ||
            std::abs(candidate.residual) < std::abs(best.residual)) {
          best = candidate;
        }
      };

      for (int k = 1; k <= 200 && !bracket_found; ++k) {
        double T_low_probe = T_guess - k * step;
        if (T_low_probe >= 50.0) {
          EvalResult low = eval_at_T(T_low_probe);
          update_best(low);
          if (low.valid && prev_low.valid &&
              low.residual * prev_low.residual <= 0.0) {
            bracket_low = low;
            bracket_high = prev_low;
            bracket_found = true;
            break;
          }
          if (low.valid) {
            prev_low = low;
          }
        }

        double T_high_probe = T_guess + k * step;
        if (T_high_probe <= 5000.0) {
          EvalResult high = eval_at_T(T_high_probe);
          update_best(high);
          if (high.valid && prev_high.valid &&
              high.residual * prev_high.residual <= 0.0) {
            bracket_low = prev_high;
            bracket_high = high;
            bracket_found = true;
            break;
          }
          if (high.valid) {
            prev_high = high;
          }
        }
      }

      if (bracket_found) {
        for (int k = 0; k < 100; ++k) {
          double T_mid = 0.5 * (bracket_low.T + bracket_high.T);
          double denom = bracket_high.residual - bracket_low.residual;
          if (std::abs(denom) > 1e-300) {
            double T_secant =
                bracket_low.T -
                bracket_low.residual * (bracket_high.T - bracket_low.T) /
                    denom;
            double T_min = std::min(bracket_low.T, bracket_high.T);
            double T_max = std::max(bracket_low.T, bracket_high.T);
            if (T_secant > T_min && T_secant < T_max) {
              T_mid = T_secant;
            }
          }

          EvalResult mid = eval_at_T(T_mid);
          update_best(mid);

          if (mid.valid && rel_error(mid.residual) <= tol) {
            best = mid;
            converged = true;
            break;
          }

          if (!mid.valid) {
            break;
          }

          if (bracket_low.residual * mid.residual <= 0.0) {
            bracket_high = mid;
          } else {
            bracket_low = mid;
          }
        }
      }

      if (best.valid) {
        T = best.T;
        h = best.h;
        if (rel_error(best.residual) <= tol) {
          converged = true;
        }
      }
    }
  }

  is_converged = converged;

  if (!converged) {
    std::cerr << "Warning: solveTemperatureFromH (VLE) did not converge "
                 "fully. Error : "
              << (h - h_target) << " Tguess " << T_guess << " Tfinal " << T
              << std::endl;
  }

  setTP(T, p);
  return T;
}

double Mixture::solvePropertiesFromH_VLE(
    double h_target, double T_guess, double tol, int max_iter,
    bool &is_converged, double &Z, double &rho, double &h, double &cp,
    double &mu, double &lambda, double &alpha, double &dpdrho_T,
    double &sound_speed, double &beta_mass) {
  bool temperature_converged = false;
  const double p = this->P;
  const std::string mixCR = "CR1";
  const double kij = 0.0;

  double T_solution = solveTemperatureFromH_VLE(
      h_target, T_guess, tol, max_iter, temperature_converged);

  setTP(T_solution, p);

  std::vector<double> x, y;
  double beta_mol = 0.0;
  std::vector<double> dummy_K;
  bool flash_success = flashTPN(this->X, this->T, this->P, mixCR,
                                this->speciesList, this->EoS, kij, x, y,
                                beta_mol, dummy_K);

  calculateProperties_VLE_fromState(
      x, y, beta_mol, flash_success, true, h_target, Z, rho, h, cp, mu,
      lambda, alpha, dpdrho_T, sound_speed, beta_mass);

  double rel_error = std::abs((h - h_target) / (h_target + 1e-10));
  is_converged = temperature_converged || rel_error <= tol;
  return T_solution;
}


// Metodo per impostare la composizione molare delle specie
void Mixture::setY(const std::vector<double> &y, const std::string &specie) {
  this->Y = y;
  mass2mol(y, this->speciesList, this->X, this->Mavg);
  this->R = R0 * 1000 / this->Mavg;

  // if (specie.empty()) {
  //     // Se non è specificata una specie, imposta tutta la composizione
  //     this->Y = y;
  //     this->X = mass2mol(y, specie)
  // } else {
  //     throw std::invalid_argument("To be implemented");
  //     // Altrimenti, aggiorna la composizione di una specie specifica
  //     // auto it = std::find(Species_names.begin(), Species_names.end(),
  //     specie);
  //     // if (it != Species_names.end()) {
  //     //     size_t index = std::distance(Species_names.begin(), it);
  //     //     this->Y[index] = y[index];
  //     // } else {
  //     //     std::cerr << "Specie " << specie << " non trovata!" <<
  //     std::endl;
  //     // }
  // }
}

// Metodo per impostare la temperatura (T) e la pressione (P)
void Mixture::setTP(double T_val, double P_val) {
  this->T = T_val;
  this->P = P_val;
}

void Mixture::setT(double T_val) { this->T = T_val; }
void Mixture::setP(double P_val) { this->P = P_val; }

void Mixture::VdWMixing(double &am, double &bm, double &dadTm, double &ddadTm) {
  VdWmixing_CR1(this->X, this->T, this->speciesList, this->EoS, 0, am, bm,
                dadTm, ddadTm); // kij = 0
}

void Mixture::VdWMixing_ambm(double &am, double &bm) {
  VdWmixing_CR1_ambm(this->X, this->T, this->speciesList, this->EoS, 0, am,
                     bm); // kij = 0
}

double Mixture::calculateZ(double am, double bm) {
  return getZ(am, bm, this->P, this->T, this->X, this->EoS, this->speciesList);
}

void Mixture::calculateZ(double &Z) {
  double am = 0.0, bm = 0.0, dadTm = 0.0, ddadTm = 0.0;
  VdWMixing(am, bm, dadTm, ddadTm);
  Z = calculateZ(am, bm);
}

double Mixture::calculateH(double am, double bm, double dadTm, double ddadTm) {
  return getH(this->T,    // temperatura
              this->P,    // pressione
              this->Mavg, // massa molare media
              this->X,    // composizione (vector<double>)
              this->Y,
              this->speciesList, // database specie
              am, bm, dadTm, ddadTm,
              this->EoS // equazione di stato
  );
} //.....

// double Mixture::calculateCp(double am, double bm, double dadTm, double
// ddadTm)
// {

// }

double Mixture::getHf() {
  // getIdealDepartureProperties(this->T, this->Y, this->Species, this->Mavg,
  // cp0, cv0, h0);
  double h, _;
  getIdealDepartureProperties(298.15, this->Y, this->speciesList, this->Mavg, h,
                              _, _);
  return h;
}

// void Mixture::calculateProperties(double& Z, double& rho, double& h, double&
// cp, double& mu, double& lambda, double& alpha)
void Mixture::calculateProperties(double &Z, double &rho, double &h, double &cp,
                                  double &mu, double &lambda, double &alpha,
                                  double &dpdrho_T, double &sound_speed)

{
  double am = 0.0, bm = 0.0, dadTm = 0.0, ddadTm = 0.0;
  VdWMixing(am, bm, dadTm, ddadTm);
  h = calculateH(am, bm, dadTm, ddadTm);
  Z = calculateZ(am, bm);

  double v;
  getRho(Z, this->T, this->P, this->Mavg, rho, v);

  double cp0, cv0, h0;
  getIdealDepartureProperties(this->T, this->Y, this->speciesList, this->Mavg,
                              cp0, cv0, h0);
  double de, dh, dcv;
  departureFunctions(am, bm, dadTm, ddadTm, this->Mavg, v, this->T, this->P,
                     this->Y, this->EoS, this->speciesList, de, dh, dcv);

  h = h0 + dh;
  double cv = cv0 + dcv;

  double dpdT_v;
  double dpdv_T;
  double dpdT_rho; // double dpdrho_T;
  thermo_derivatives(am, bm, dadTm, v, rho, this->Mavg, this->T, this->X,
                     this->speciesList, dpdT_v, dpdv_T, dpdT_rho, dpdrho_T);
  cp = cv + (T / (rho * rho)) * ((dpdT_rho * dpdT_rho) / dpdrho_T);

  sound_speed =
      std::sqrt(dpdrho_T + (this->T / (cv)) * std::pow(dpdT_rho / rho, 2));

  // double Tc_m; double Vc_m; double omega_m; double eps_m; double M_m; double
  // Fcm;
  double sigma3m, epsNm, omegaNm, MNm, Vc_m, eps_m, omega_m, M_m, Fcm, Tc_m;

  chung_mixing(this->X, this->speciesList, sigma3m, epsNm, omegaNm, MNm, Vc_m,
               eps_m, omega_m, M_m, Fcm, Tc_m);

  if (omega_m < 0) {
    omega_m = 0;
  }

  if (Fcm > 1) {
    Fcm = 1;
  }
  Fcm = 1;

  if (cp < 1e-4)
    cp = 1e-4; // Stability floor

  double mu0;
  chung_transport(eps_m, Fcm, M_m, Vc_m, omega_m, Tc_m, cv0, v, this->T, mu,
                  lambda, mu0);
  alpha = lambda / cp;
  if (alpha < 0)
    alpha = 1e-8; // Ensure positive diffusivity
}

void Mixture::calculatePropertiesForPhase(const std::string &phase, double &Z,
                                          double &rho, double &h, double &cp,
                                          double &mu, double &lambda,
                                          double &alpha, double &dpdrho_T,
                                          double &sound_speed) {
  if (phase != "Liq" && phase != "Vap") {
    throw std::invalid_argument("phase must be either 'Liq' or 'Vap'");
  }

  double am = 0.0, bm = 0.0, dadTm = 0.0, ddadTm = 0.0;
  VdWMixing(am, bm, dadTm, ddadTm);
  Z = getZ(am, bm, this->P, this->T, this->X, phase, this->speciesList);

  double v;
  getRho(Z, this->T, this->P, this->Mavg, rho, v);

  double cp0, cv0, h0;
  getIdealDepartureProperties(this->T, this->Y, this->speciesList, this->Mavg,
                              cp0, cv0, h0);
  double de, dh, dcv;
  departureFunctions(am, bm, dadTm, ddadTm, this->Mavg, v, this->T, this->P,
                     this->Y, this->EoS, this->speciesList, de, dh, dcv);

  h = h0 + dh;
  double cv = cv0 + dcv;

  double dpdT_v;
  double dpdv_T;
  double dpdT_rho;
  thermo_derivatives(am, bm, dadTm, v, rho, this->Mavg, this->T, this->X,
                     this->speciesList, dpdT_v, dpdv_T, dpdT_rho, dpdrho_T);
  cp = cv + (T / (rho * rho)) * ((dpdT_rho * dpdT_rho) / dpdrho_T);

  sound_speed =
      std::sqrt(dpdrho_T + (this->T / (cv)) * std::pow(dpdT_rho / rho, 2));

  double sigma3m, epsNm, omegaNm, MNm, Vc_m, eps_m, omega_m, M_m, Fcm, Tc_m;

  chung_mixing(this->X, this->speciesList, sigma3m, epsNm, omegaNm, MNm, Vc_m,
               eps_m, omega_m, M_m, Fcm, Tc_m);

  if (omega_m < 0) {
    omega_m = 0;
  }

  if (Fcm > 1) {
    Fcm = 1;
  }
  Fcm = 1;

  if (cp < 1e-4)
    cp = 1e-4;

  double mu0;
  chung_transport(eps_m, Fcm, M_m, Vc_m, omega_m, Tc_m, cv0, v, this->T, mu,
                  lambda, mu0);
  alpha = lambda / cp;
  if (alpha < 0)
    alpha = 1e-8;
}

double Mixture::calculateCriticalTemperature() {
  // sigma3m, epsNm, omegaNm, MNm, Vc_m, eps_m, omega_m, M_m, Fcm, Tc_m
  double sigma3m, epsNm, omegaNm, MNm, Vc_m, eps_m, omega_m, M_m, Fcm, Tc_m;
  chung_mixing(this->X, this->speciesList, sigma3m, epsNm, omegaNm, MNm, Vc_m,
               eps_m, omega_m, M_m, Fcm, Tc_m);
  return Tc_m;
}

double Mixture::calculateCriticalPressure() {
  double sigma3m, epsNm, omegaNm, MNm, Vc_m, eps_m, omega_m, M_m, Fcm, Tc_m;
  chung_mixing(this->X, this->speciesList, sigma3m, epsNm, omegaNm, MNm, Vc_m,
               eps_m, omega_m, M_m, Fcm, Tc_m);

  double Zcm = 0.0;
  for (size_t i = 0; i < this->X.size(); ++i) {
    Zcm += this->X[i] * this->speciesList[i].Zc;
  }

  if (Vc_m <= 0.0)
    return 0.0;

  // Pcm = Zcm * R0 * Tcm / Vcm
  // R0 is 8.314472 J/(mol*K)
  // Tc_m in K, Vc_m in m^3/mol
  return (Zcm * 8.314472 * Tc_m) / Vc_m;
}

// double getH(
//     double T, double p, double M,
//     const std::vector<double>& Y,
//     const std::vector<Species>& species,
//     double Z,
//     double am, double bm, double dadTm, double ddadTm,
//     const std::string& EoS
// ) {
//     double cp0, cv0, h0;
//     getIdealDepartureProperties(T, Y, species, M, cp0, cv0, h0);

//     double v;
//     double rho;
//     getRho(Z, T, p, M, rho, v);

//     double de, dh, dcv;
//     departureFunctions(am, bm, dadTm, ddadTm, M, v, T, p, Y, EoS, species,
//     de, dh, dcv);

//     return h0 + dh;
// }

double Mixture::calculateRho(double am, double bm) {
  return this->P / (getZ(am, bm, this->P, this->T, this->X, this->EoS,
                         this->speciesList) *
                    this->R * this->T);
}

double Mixture::calculateRho(double Z) {
  return this->P / (Z * this->R * this->T);
}

bool Mixture::checkStability(double &TPDmin, std::vector<double> &K_init,
                             const std::string &mixCR, double kij) {
  bool stable = true;
  TPD_SSI(this->X, this->T, this->P, mixCR, this->speciesList, this->EoS, kij,
          stable, TPDmin, K_init);
  return stable;
}

bool Mixture::calculateBeta(double &beta_mass, double &beta_mol,
                            std::vector<double> &x, std::vector<double> &y,
                            const std::string &mixCR, double kij) {
  std::vector<double> dummy_K_init; // No initial K factors provided directly
  bool success = flashTPN(this->X, this->T, this->P, mixCR, this->speciesList,
                          this->EoS, kij, x, y, beta_mol, dummy_K_init);

  if (success) {
    double M_l = 0.0;
    double M_v = 0.0;
    for (size_t i = 0; i < this->X.size(); ++i) {
      M_l += x[i] * this->speciesList[i].M;
      M_v += y[i] * this->speciesList[i].M;
    }

    double M_mix = 0.0;
    for (size_t i = 0; i < this->X.size(); ++i) {
      M_mix += this->X[i] * this->speciesList[i].M;
    } // Overall mixture molar mass

    beta_mass = beta_mol * (M_v / M_mix);
  } else {
    beta_mass = 0.0;
  }

  return success;
}

bool Mixture::checkSpinodal(double &detQ, const std::string &mixCR,
                            double kij) {
  return ::checkSpinodal(this->X, this->T, this->P, mixCR, this->speciesList,
                         this->EoS, kij, detQ);
}

// Calcola le proprietà termodinamiche e di trasporto in condizioni VLE.
// Segue l'algoritmo RFVLEproperties.m descritto in report_VLE.md.
void Mixture::calculateProperties_VLE(double &Z, double &rho, double &h,
                                      double &cp, double &mu, double &lambda,
                                      double &alpha, double &dpdrho_T,
                                      double &sound_speed, double &beta_mass) {
  const std::string mixCR = "CR1";
  const double kij = 0.0;

  std::vector<double> x, y;
  double beta_mol = 0.0;
  std::vector<double> dummy_K;
  bool success = flashTPN(this->X, this->T, this->P, mixCR, this->speciesList,
                          this->EoS, kij, x, y, beta_mol, dummy_K);

  calculateProperties_VLE_fromState(
      x, y, beta_mol, success, false, 0.0, Z, rho, h, cp, mu, lambda, alpha,
      dpdrho_T, sound_speed, beta_mass);
}

void Mixture::calculateProperties_VLE_fromState(
    const std::vector<double> &x, const std::vector<double> &y,
    double beta_mol, bool flash_success, bool force_enthalpy_beta,
    double h_target, double &Z, double &rho, double &h, double &cp, double &mu,
    double &lambda, double &alpha, double &dpdrho_T, double &sound_speed,
    double &beta_mass) {
  const double kij = 0.0;

  // --- 2) Masse molari delle fasi ---
  double M_l = 0.0, M_v = 0.0;
  for (size_t i = 0; i < this->X.size(); ++i) {
    M_l += x[i] * this->speciesList[i].M;
    M_v += y[i] * this->speciesList[i].M;
  }

  // Frazione massica di vapore: beta_mass = beta_mol * (M_v / M_mix)
  beta_mass = (this->Mavg > 0.0) ? beta_mol * (M_v / this->Mavg) : 0.0;

  // --- 3) Caso monofase: delega a calculateProperties ---
  // Se beta_mass è 0 (tutto liquido) o 1 (tutto vapore) non c'è bifase.
  if (!flash_success ||
      (!force_enthalpy_beta && (beta_mass <= 0.0 || beta_mass >= 1.0))) {
    if (beta_mass < 0.0) beta_mass = 0.0;
    if (beta_mass > 1.0) beta_mass = 1.0;
    calculateProperties(Z, rho, h, cp, mu, lambda, alpha, dpdrho_T,
                        sound_speed);
    return;
  }

  // --- 4) Composizioni massiche per fase ---
  std::vector<double> x_mass(this->X.size(), 0.0);
  std::vector<double> y_mass(this->X.size(), 0.0);
  for (size_t i = 0; i < this->X.size(); ++i) {
    if (M_l > 0.0) x_mass[i] = (x[i] * this->speciesList[i].M) / M_l;
    if (M_v > 0.0) y_mass[i] = (y[i] * this->speciesList[i].M) / M_v;
  }

  // --- 5) Regole di miscelazione VdW per le due fasi ---
  double am_l, bm_l, dadTm_l, ddadTm_l;
  double am_v, bm_v, dadTm_v, ddadTm_v;
  VdWmixing_CR1(x, this->T, this->speciesList, this->EoS, kij, am_l, bm_l,
                dadTm_l, ddadTm_l);
  VdWmixing_CR1(y, this->T, this->speciesList, this->EoS, kij, am_v, bm_v,
                dadTm_v, ddadTm_v);

  // --- 6) Fattori di compressibilità e volumi molari ---
  double Z_l = getZ(am_l, bm_l, this->P, this->T, x, "Liq", this->speciesList);
  double Z_v = getZ(am_v, bm_v, this->P, this->T, y, "Vap", this->speciesList);

  double v_l = (Z_l * R0 * this->T) / this->P; // m³/mol
  double v_v = (Z_v * R0 * this->T) / this->P;

  double rho_l = M_l / v_l; // kg/m³
  double rho_v = M_v / v_v;

  // --- 7) Proprietà ideali per fase (NASA7 via getIdealDepartureProperties) ---
  double cp0_l, cv0_l, h0_l;
  double cp0_v, cv0_v, h0_v;
  getIdealDepartureProperties(this->T, x_mass, this->speciesList, M_l, cp0_l,
                              cv0_l, h0_l);
  getIdealDepartureProperties(this->T, y_mass, this->speciesList, M_v, cp0_v,
                              cv0_v, h0_v);

  // cv0 è già calcolato correttamente da getIdealDepartureProperties
  // come cp0 - R0/M (con R0 e M in unità coerenti). NON sovrascrivere.

  // --- 8) Funzioni di scostamento ---
  double de_l, dh_l = 0.0, dcv_l = 0.0;
  double de_v, dh_v = 0.0, dcv_v = 0.0;
  departureFunctions(am_l, bm_l, dadTm_l, ddadTm_l, M_l, v_l, this->T,
                     this->P, x_mass, this->EoS, this->speciesList, de_l, dh_l,
                     dcv_l);
  departureFunctions(am_v, bm_v, dadTm_v, ddadTm_v, M_v, v_v, this->T,
                     this->P, y_mass, this->EoS, this->speciesList, de_v, dh_v,
                     dcv_v);

  // --- 9) Entalpia e cv reali per fase ---
  double h_l = h0_l + dh_l;
  double h_v = h0_v + dh_v;
  double cv_l = cv0_l + dcv_l;
  double cv_v = cv0_v + dcv_v;

  if (force_enthalpy_beta && std::isfinite(h_target) &&
      std::abs(h_v - h_l) > 1e-12) {
    double beta_from_h = (h_target - h_l) / (h_v - h_l);
    if (beta_from_h >= 0.0 && beta_from_h <= 1.0) {
      beta_mass = beta_from_h;
    }
  }

  if (beta_mass <= 0.0 || beta_mass >= 1.0) {
    if (beta_mass < 0.0) beta_mass = 0.0;
    if (beta_mass > 1.0) beta_mass = 1.0;
    calculateProperties(Z, rho, h, cp, mu, lambda, alpha, dpdrho_T,
                        sound_speed);
    return;
  }

  // --- 10) Derivate termodinamiche per Cp ---
  double dpdT_v_l, dpdv_T_l, dpdT_rho_l, dpdrho_T_l;
  double dpdT_v_v, dpdv_T_v, dpdT_rho_v, dpdrho_T_v;
  thermo_derivatives(am_l, bm_l, dadTm_l, v_l, rho_l, M_l, this->T, x,
                     this->speciesList, dpdT_v_l, dpdv_T_l, dpdT_rho_l,
                     dpdrho_T_l);
  thermo_derivatives(am_v, bm_v, dadTm_v, v_v, rho_v, M_v, this->T, y,
                     this->speciesList, dpdT_v_v, dpdv_T_v, dpdT_rho_v,
                     dpdrho_T_v);

  // cp = cv - T*(dp/dT)_v^2 / (dp/dv)_T / M
  // (dp/dv)_T < 0 per fluido stabile, quindi cp > cv come atteso.
  double cp_l = (dpdv_T_l != 0.0)
                    ? cv_l - (this->T * dpdT_v_l * dpdT_v_l / dpdv_T_l) / M_l
                    : cv_l;
  double cp_v = (dpdv_T_v != 0.0)
                    ? cv_v - (this->T * dpdT_v_v * dpdT_v_v / dpdv_T_v) / M_v
                    : cv_v;

  if (cp_l < 1e-4) cp_l = 1e-4;
  if (cp_v < 1e-4) cp_v = 1e-4;

  // --- 11) Velocità del suono per fase ---
  // SEMPLIFICAZIONE: formula monofase applicata a ciascuna fase separatamente.
  // La velocità del suono bifase rigorosa richiederebbe la formula di Wood,
  // che dipende dalla comprimibilità del sistema bifase completo.
  double c_l = std::sqrt(std::max(0.0, dpdrho_T_l +
                         (this->T / cv_l) * std::pow(dpdT_rho_l / rho_l, 2)));
  double c_v = std::sqrt(std::max(0.0, dpdrho_T_v +
                         (this->T / cv_v) * std::pow(dpdT_rho_v / rho_v, 2)));

  // --- 12) Proprietà di trasporto per fase (Chung) ---
  // SEMPLIFICAZIONE: le regole di Chung vengono applicate usando le frazioni
  // molari di fase (x, y) invece della composizione globale. Le proprietà
  // critiche di miscela sono quindi specifiche di ogni fase.
  double sigma3m, epsNm, omegaNm, MNm, Vc_m, eps_m, omega_m, M_m, Fcm, Tc_m;

  chung_mixing(x, this->speciesList, sigma3m, epsNm, omegaNm, MNm, Vc_m,
               eps_m, omega_m, M_m, Fcm, Tc_m);
  if (omega_m < 0.0) omega_m = 0.0;
  Fcm = 1.0; // coerente con calculateProperties
  double mu_l, lambda_l, mu0_l;
  chung_transport(eps_m, Fcm, M_m, Vc_m, omega_m, Tc_m, cv0_l, v_l, this->T,
                  mu_l, lambda_l, mu0_l);

  chung_mixing(y, this->speciesList, sigma3m, epsNm, omegaNm, MNm, Vc_m,
               eps_m, omega_m, M_m, Fcm, Tc_m);
  if (omega_m < 0.0) omega_m = 0.0;
  Fcm = 1.0;
  double mu_v, lambda_v, mu0_v;
  chung_transport(eps_m, Fcm, M_m, Vc_m, omega_m, Tc_m, cv0_v, v_v, this->T,
                  mu_v, lambda_v, mu0_v);

  // --- 13) Medie massiche per la miscela bifase ---
  h  = beta_mass * h_v  + (1.0 - beta_mass) * h_l;
  cp = beta_mass * cp_v + (1.0 - beta_mass) * cp_l;
  mu = beta_mass * mu_v + (1.0 - beta_mass) * mu_l;
  double lambda_mix = beta_mass * lambda_v + (1.0 - beta_mass) * lambda_l;
  lambda = lambda_mix;

  // SEMPLIFICAZIONE: rho bifase dalla media volumetrica (regola delle lever).
  // 1/rho = beta_mass/rho_v + (1-beta_mass)/rho_l
  double inv_rho = beta_mass / rho_v + (1.0 - beta_mass) / rho_l;
  rho = 1.0 / inv_rho;

  // Z dalla densità risultante: Z = P / (rho * R_spec * T)
  // SEMPLIFICAZIONE: Z effettivo bifase, non ha significato fisico diretto
  // come in monofase ma è necessario per l'interfaccia.
  double R_spec = R0 * 1000.0 / this->Mavg; // J/(kg·K)
  Z = this->P / (rho * R_spec * this->T);

  // dpdrho_T bifase: media massica delle derivate delle singole fasi.
  // SEMPLIFICAZIONE: la derivata rigorosa bifase richiederebbe la Jacobiana
  // completa del sistema flash.
  dpdrho_T = beta_mass * dpdrho_T_v + (1.0 - beta_mass) * dpdrho_T_l;

  // Velocità del suono bifase: media massica.
  // SEMPLIFICAZIONE: la formula di Wood sarebbe
  //   1/(rho*c²) = beta_mass/(rho_v*c_v²) + (1-beta_mass)/(rho_l*c_l²)
  // ma è più costosa e richiede c_l, c_v > 0 garantiti.
  sound_speed = beta_mass * c_v + (1.0 - beta_mass) * c_l;

  alpha = lambda / cp;
  if (alpha < 0.0) alpha = 1e-8;
}
