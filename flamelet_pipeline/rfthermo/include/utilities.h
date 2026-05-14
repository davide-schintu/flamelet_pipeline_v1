#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <string>
#include <tuple>

struct Species
{
    std::string Names;
    double Tc;
    double Pc;
    double omega;
    double M;
    double Vc;
    double Zc;
    double vd;
    double AlowT, BlowT, ClowT, DlowT, ElowT, FlowT, GlowT;
    double AHigT, BHigT, CHigT, DHigT, EHigT, FHigT, GHigT;
    double Tlow, Thig, Tmed;
    double geometry;
    double eps_kb;
    double sigma;
    double mu;
    double alpha;
    double Zrot;
};

void mass2mol(const std::vector<double>& Y, const std::vector<Species>& Species, std::vector<double>& X, double& Mw);
std::vector<Species> database_reader(const std::string& data_dir, const std::vector<std::string>& orderedSpeciesNames = {});



#endif // UTILITIES_H
