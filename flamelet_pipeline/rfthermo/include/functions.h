#ifndef EOS_FUNCTIONS_H
#define EOS_FUNCTIONS_H

#include "mixingRules.h"
#include <map>
#include <string>
#include <vector>

double getZ(double am, double bm, double P, double T,
            const std::vector<double> &X, const std::string &EoS,
            const std::vector<Species> &species);

double getH(double T, double p, double M, const std::vector<double> &Y,
            const std::vector<double> &X, const std::vector<Species> &species,
            double am, double bm, double dadTm, double ddadTm,
            const std::string &EoS);

void getRho(double Z, double T, double p, double M, double &rho, double &v);

void getIdealDepartureProperties(double T, const std::vector<double> &Y,
                                 const std::vector<Species> &species, double M,
                                 double &cp0, double &cv0, double &h0);

void departureFunctions(double am, double bm, double dadTm, double ddadTm,
                        double M, double v, double T, double p,
                        const std::vector<double> &X, const std::string &EoS,
                        const std::vector<Species> &species,
                        double &de, // output
                        double &dh, // output
                        double &dcv // output
);

void thermo_derivatives(double am, double bm, double dadTm, double v,
                        double rho, double M, double T,
                        const std::vector<double> &X,
                        const std::vector<Species> &species, double &dpdT_v,
                        double &dpdv_T, double &dpdT_rho, double &dpdrho_T);

double collision_integrale(double);

void chung_mixing(const std::vector<double> &X,
                  const std::vector<Species> &Species, double &sigma3m,
                  double &epsNm, double &omegaNm, double &MNm, double &Vc_m,
                  double &eps_m, double &omega_m, double &M_m, double &Fcm,
                  double &Tc_m);

void chung_transport(double eps_m, double Fcm, double W_m, double Vc_m,
                     double omega_m, double Tc_m, double cv0, double v,
                     double T, double &mu, double &lambda_, double &mu0);

double fug_coef(const std::vector<double> &X, double bi, double bm,
                double sum_ref_i, double am, double Z, double T, double p,
                const std::string &EoS, const std::vector<Species> &species);

void TPD_SSI(const std::vector<double> &X, double T, double p,
             const std::string &mixCR, const std::vector<Species> &Species,
             const std::string &EoS, double kij, bool &stable, double &TPDmin,
             std::vector<double> &K_init_min);

bool flashTPN(const std::vector<double> &z, double T, double p,
              const std::string &mixCR, const std::vector<Species> &SpeciesList,
              const std::string &EoS, double kij, std::vector<double> &x,
              std::vector<double> &y, double &psi_v,
              const std::vector<double> &K_init);

// Spinodal stability check via Hessian of Gibbs energy
// Returns: true if mixture is inside spinodal (diffusionally unstable)
//          false if outside spinodal (metastable or stable)
// detQ:   the determinant of the stability matrix Q (output)
bool checkSpinodal(const std::vector<double> &X, double T, double p,
                   const std::string &mixCR,
                   const std::vector<Species> &Species, const std::string &EoS,
                   double kij, double &detQ);

#endif // EOS_FUNCTIONS_H
