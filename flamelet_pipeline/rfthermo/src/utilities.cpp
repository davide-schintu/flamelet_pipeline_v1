#include "../include/utilities.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <iomanip> // For formatting output
#include <string> // Explicitly include string

const double R0 = 8.314472;

void mass2mol(const std::vector<double>& Y, const std::vector<Species>& Species, std::vector<double>& X, double& Mavg) {
    int num_species = Y.size();
    std::vector<double> oneM(num_species);
    
    for (int i = 0; i < num_species; ++i) {
        oneM[i] = Y[i] / Species[i].M;
    }

    double sum_oneM = std::accumulate(oneM.begin(), oneM.end(), 0.0);
    Mavg = (sum_oneM == 0.0) ? 0.0 : 1.0 / sum_oneM;

    for (int i = 0; i < num_species; ++i) {
        X[i] = oneM[i] * Mavg;
    }

    double sum_X = std::accumulate(X.begin(), X.end(), 0.0);
    if (Mavg != 0.0 && std::abs(sum_X - 1.0) > 1e-9) {
        for (int i = 0; i < num_species; ++i) {
            X[i] /= sum_X;
        }
    }
}

#include <map>
#include <iostream>
#include <algorithm>

std::vector<Species> database_reader(const std::string& filename, const std::vector<std::string>& orderedSpeciesNames)
{
    std::map<std::string, Species> species_map;
    std::vector<Species> species_list;
    std::ifstream file(filename);
    
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return species_list;
    }

    std::string line;
    std::getline(file, line); // Skip header

    const double R0 = 8.314472;

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        if (tokens.size() < 31) {
            std::cerr << "Linea malformata: " << line << std::endl;
            continue;
        }

        Species s;
        size_t i = 0;

        s.Names = tokens[i++];
        s.Tc = std::stod(tokens[i++]);
        s.Pc = std::stod(tokens[i++]);
        s.omega = std::stod(tokens[i++]);
        s.M = std::stod(tokens[i++]);
        s.Vc = std::stod(tokens[i++]);
        s.Zc = std::stod(tokens[i++]);
        s.vd = std::stod(tokens[i++]);

        s.AlowT = std::stod(tokens[i++]);
        s.BlowT = std::stod(tokens[i++]);
        s.ClowT = std::stod(tokens[i++]);
        s.DlowT = std::stod(tokens[i++]);
        s.ElowT = std::stod(tokens[i++]);
        s.FlowT = std::stod(tokens[i++]);
        s.GlowT = std::stod(tokens[i++]);

        s.AHigT = std::stod(tokens[i++]);
        s.BHigT = std::stod(tokens[i++]);
        s.CHigT = std::stod(tokens[i++]);
        s.DHigT = std::stod(tokens[i++]);
        s.EHigT = std::stod(tokens[i++]);
        s.FHigT = std::stod(tokens[i++]);
        s.GHigT = std::stod(tokens[i++]);

        s.Tlow = std::stod(tokens[i++]);
        s.Thig = std::stod(tokens[i++]);
        s.Tmed = std::stod(tokens[i++]);

        s.geometry = std::stoi(tokens[i++]);
        s.eps_kb = std::stod(tokens[i++]);
        s.sigma = std::stod(tokens[i++]);
        s.mu = std::stod(tokens[i++]);
        s.alpha = std::stod(tokens[i++]);
        s.Zrot = std::stod(tokens[i++]);

        // Conversioni (se necessarie)
        // Conversioni (se necessarie)
        s.Zc = (s.Tc != 0.0) ? (s.Pc * s.Vc) / (R0 * s.Tc) : 1.0;

        // If no ordering required, just push back
        if (orderedSpeciesNames.empty()) {
            species_list.push_back(s);
        } else {
             // Store in map for reordering later
             species_map[s.Names] = s;
        }
    }

    // If ordering required, rebuild the list in order
    if (!orderedSpeciesNames.empty()) {
        for (const auto& name : orderedSpeciesNames) {
            auto it = species_map.find(name);
            if (it != species_map.end()) {
                species_list.push_back(it->second);
            } else {
                std::cerr << "CRITICAL ERROR: Species '" << name << "' requested by OpenFOAM but NOT FOUND in " << filename << std::endl;
                // Optionally throw or exit, but for now strict warning + empty species might crash
                exit(1);
            }
        }
    }

    file.close();
    return species_list;
}
