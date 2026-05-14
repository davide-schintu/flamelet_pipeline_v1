#ifndef MIXTUREPROPERTIES_H
#define MIXTUREPROPERTIES_H

#include <vector>
#include <string>

// Serve la definizione di Species
#include "utilities.h"

// Costante universale dei gas
constexpr double R0 = 8.314472;

// Dichiarazioni delle funzioni
void SRK_abc(double Tc, double Pc, double omega, double eps_kb, double sigma,
             double& ac, double& bc, double& cc);

void alphaFunction(double T, double Tc, double cc, const std::string& EoS,
                   double& alpha, double& dalphadT, double& ddalphadT);

void VdWmixing_CR1(const std::vector<double>& X, double T,
                   const std::vector<Species>& Species, const std::string& EoS, double kij,
                   double& am, double& bm, double& dadTm, double& ddadTm);

void VdWmixing_CR1_ambm(const std::vector<double>& X, double T, 
                   const std::vector<Species>& Species, const std::string& EoS, double kij, 
                   double& am, double& bm);

void VdWmixing_CR2(const std::vector<double>& X, double T,
                   const std::vector<Species>& Species, const std::string& EoS, double kij,
                   double& am, double& bm, double& dadTm, double& ddadTm);

void VdWmixing_CR2_ambm(const std::vector<double>& X, double T, 
                   const std::vector<Species>& Species, const std::string& EoS, double kij, 
                   double& am, double& bm);

void VdWmixing_CR1_fugacity(const std::vector<double>& X, double T,
                            const std::vector<Species>& Species, const std::string& EoS, double kij,
                            double& am, double& bm, double& dadTm, double& ddadTm,
                            std::vector<double>& sum_ref_i, std::vector<double>& bc);

void VdWmixing_CR2_fugacity(const std::vector<double>& X, double T,
                            const std::vector<Species>& Species, const std::string& EoS, double kij,
                            double& am, double& bm, double& dadTm, double& ddadTm,
                            std::vector<double>& sum_ref_i, std::vector<double>& bc);

#endif // MIXTUREPROPERTIES_H
