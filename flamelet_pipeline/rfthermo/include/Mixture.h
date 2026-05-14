#ifndef MIXTURE_H
#define MIXTURE_H

#include "functions.h"
#include <iostream>
#include <string>
#include <vector>

// Dichiarazione della classe Mixture
class Mixture {
public:
  // Costruttore che inizializza il database
  Mixture(std::string database_name,
          const std::vector<std::string> &orderedSpeciesNames = {});

  // Metodo per impostare la composizione Y
  void setY(const std::vector<double> &y, const std::string &specie = "");

  void setT(double T_val);
  void setP(double P_val);
  void setTP(double, double);

  void VdWMixing(double &am, double &bm, double &dadTm, double &ddadTm);
  void VdWMixing_ambm(double &am, double &bm);

  // Get properties. Performances must be verified
  double calculateZ(double am, double bm);
  void calculateZ(double &Z);
  double calculateH(double am, double bm, double dadTm, double ddadTm);
  double getHf();

  void calculateProperties(double &Z, double &rho, double &h, double &cp,
                           double &mu, double &lambda, double &alpha,
                           double &dpdrho_T, double &sound_speed);
  void calculatePropertiesForPhase(const std::string &phase, double &Z,
                                   double &rho, double &h, double &cp,
                                   double &mu, double &lambda, double &alpha,
                                   double &dpdrho_T, double &sound_speed);

  double calculateCriticalTemperature();
  double calculateCriticalPressure();

  // double solveTemperatureFromH(double h_target,
  //                          double T_low, double T_high, double tol, int
  //                          max_iter);

  double solveTemperatureFromH(double h_target, double T_guess, double tol,
                               int max_iter);

  double solveTemperatureFromH_VLE(double h_target, double T_guess, double tol,
                               int max_iter, bool &converged);

  double calculateRho(double Z);
  double calculateRho(double am, double bm);

  bool checkStability(double &TPDmin, std::vector<double> &K_init,
                      const std::string &mixCR, double kij);
  bool calculateBeta(double &beta_mass, double &beta_mol,
                     std::vector<double> &x, std::vector<double> &y,
                     const std::string &mixCR, double kij);

  bool checkSpinodal(double &detQ, const std::string &mixCR, double kij);

  void calculateProperties_VLE(double &Z, double &rho, double &h, double &cp,
                               double &mu, double &lambda, double &alpha,
                               double &dpdrho_T, double &sound_speed,
                               double &beta_mass);
  void calculateProperties_VLE_fromState(
      const std::vector<double> &x, const std::vector<double> &y,
      double beta_mol, bool flash_success, bool force_enthalpy_beta,
      double h_target, double &Z, double &rho, double &h, double &cp,
      double &mu, double &lambda, double &alpha, double &dpdrho_T,
      double &sound_speed, double &beta_mass);
  double solvePropertiesFromH_VLE(double h_target, double T_guess, double tol,
                                  int max_iter, bool &converged, double &Z,
                                  double &rho, double &h, double &cp,
                                  double &mu, double &lambda, double &alpha,
                                  double &dpdrho_T, double &sound_speed,
                                  double &beta_mass);
  double getR() const { return R; }

  //... ecc.

  // Variabili membro (attributi della classe)
  std::vector<Species> speciesList;       // Lista delle specie
  std::vector<std::string> Species_names; // Nomi delle specie
  std::vector<double> Y;                  // Mass fraction
  std::vector<double> X;                  // Mole fraction

  double T; // Temperatura
  double P; // Pressione
  double Mavg;
  double R;
  std::string EoS;
};

#endif // MIXTURE_H
