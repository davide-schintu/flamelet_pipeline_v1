// #include "../include/mixingRules.h"
// #include "../include/functions.h"
#include "../include/Mixture.h"
#include <iostream>
#include <vector>
#include <chrono>


int main() {
    
    Mixture mix = Mixture("data/ThermoProperties.csv");
    auto start = std::chrono::high_resolution_clock::now();

    int Ns = mix.Species.size();
    std::vector<double> Y(Ns, 0.0); Y[0] = 0; Y[3] = 1; // 


    mix.setY(Y);
    mix.setTP(100, 100e5);

    // void Mixture::calculateProperties(double& Z, double& rho, double& h, double& cp, double& mu, double& lambda, double& alpha)
    double Z, rho, h, cp, mu, lambda, alpha;

    mix.calculateProperties(Z, rho, h, cp, mu, lambda, alpha);

    std::cout << "Z: " << Z << std::endl;
    std::cout << "Rho: " << rho << std::endl;
    std::cout << "H: " << h << std::endl;
    std::cout << "Cp: " << cp << std::endl;
    std::cout << "Mu: " << mu << std::endl;
    std::cout << "Lambda: " << lambda << std::endl;
    std::cout << "Alpha: " << alpha << std::endl;


    std::cout << "Hf: " << mix.getHf() << std::endl;
    std::cout << "Err Hf: " << (mix.getHf() + h)/(-h) << std::endl;

    //

    // double h_target = -387996; // esempio in J/kg
    // mix.setTP(150, 80e5);

    // try {
    //     double T = solveTemperatureFromH(mix, Y, 80e5, h_target);
    //     std::cout << "Trovata temperatura: " << T << " K" << std::endl;
    // } catch (const std::exception& e) {
    //     std::cerr << e.what() << std::endl;
    // }


    return 0;
}