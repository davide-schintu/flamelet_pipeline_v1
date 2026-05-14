#include "../include/Mixture.h"
#include "../include/utilities.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(RFThermo, m) {
  m.doc() = "RFThermoFunctions Python Wrapper";
  m.attr("__build_marker__") =
      "RFThermo debug marker 2026-05-11 HP-props-v3";

  py::class_<Species>(m, "Species")
      .def_readonly("Names", &Species::Names)
      .def_readonly("Tc", &Species::Tc)
      .def_readonly("Pc", &Species::Pc)
      .def_readonly("M", &Species::M);
  // Add other fields as needed

  py::class_<Mixture>(m, "Mixture")
      .def(py::init<std::string, const std::vector<std::string> &>(),
           py::arg("database_name"),
           py::arg("orderedSpeciesNames") = std::vector<std::string>{})
      .def("setY", &Mixture::setY, py::arg("y"), py::arg("specie") = "")
      .def("setTP", &Mixture::setTP)
      .def("setT", &Mixture::setT)
      .def("setP", &Mixture::setP)
      .def("getHf", &Mixture::getHf)
      .def("calculateProperties",
           [](Mixture &self) {
             double Z, rho, h, cp, mu, lambda, alpha, dpdrho_T, sound_speed;
             self.calculateProperties(Z, rho, h, cp, mu, lambda, alpha,
                                      dpdrho_T, sound_speed);
             return py::make_tuple(Z, rho, h, cp, mu, lambda, alpha, dpdrho_T,
                                   sound_speed);
           })
      .def("calculatePropertiesForPhase",
           [](Mixture &self, const std::string &phase) {
             double Z, rho, h, cp, mu, lambda, alpha, dpdrho_T, sound_speed;
             self.calculatePropertiesForPhase(phase, Z, rho, h, cp, mu, lambda,
                                              alpha, dpdrho_T, sound_speed);
             return py::make_tuple(Z, rho, h, cp, mu, lambda, alpha, dpdrho_T,
                                   sound_speed);
           })
      .def("calculateCriticalTemperature",
           &Mixture::calculateCriticalTemperature)
      .def("calculateCriticalPressure", &Mixture::calculateCriticalPressure)

      .def("solveTemperatureFromH", &Mixture::solveTemperatureFromH)
      .def("solveTemperatureFromH_VLE",
           [](Mixture &self, double h_target, double T_guess, double tol,
              int max_iter) {
             bool is_converged = false;
             double T = self.solveTemperatureFromH_VLE(h_target, T_guess, tol,
                                                       max_iter, is_converged);
             return py::make_tuple(T, is_converged);
           },
           py::arg("h_target"), py::arg("T_guess"),
           py::arg("tol") = 1e-6, py::arg("max_iter") = 100)
      .def("calculateProperties_VLE",
           [](Mixture &self) {
             double Z, rho, h, cp, mu, lambda, alpha, dpdrho_T, sound_speed,
                 beta_mass;
             self.calculateProperties_VLE(Z, rho, h, cp, mu, lambda, alpha,
                                          dpdrho_T, sound_speed, beta_mass);
             return py::make_tuple(Z, rho, h, cp, mu, lambda, alpha, dpdrho_T,
                                   sound_speed, beta_mass);
           })
      .def("solvePropertiesFromH_VLE",
           [](Mixture &self, double h_target, double T_guess, double tol,
              int max_iter) {
             bool is_converged = false;
             double Z, rho, h, cp, mu, lambda, alpha, dpdrho_T, sound_speed,
                 beta_mass;
             double T = self.solvePropertiesFromH_VLE(
                 h_target, T_guess, tol, max_iter, is_converged, Z, rho, h, cp,
                 mu, lambda, alpha, dpdrho_T, sound_speed, beta_mass);
             return py::make_tuple(T, is_converged, Z, rho, h, cp, mu, lambda,
                                   alpha, dpdrho_T, sound_speed, beta_mass);
           },
           py::arg("h_target"), py::arg("T_guess"),
           py::arg("tol") = 1e-6, py::arg("max_iter") = 100)
      .def("calculateRho", py::overload_cast<double>(&Mixture::calculateRho))
      .def(
          "checkStability",
          [](Mixture &self, const std::string &mixCR, double kij) {
            double TPDmin = 0.0;
            std::vector<double> K_init;
            bool stable = self.checkStability(TPDmin, K_init, mixCR, kij);
            return py::make_tuple(stable, TPDmin, K_init);
          },
          py::arg("mixCR") = "CR2", py::arg("kij") = 0.0)
      .def(
          "calculateBeta",
          [](Mixture &self, const std::string &mixCR, double kij) {
            double beta_mass = 0.0;
            double beta_mol = 0.0;
            std::vector<double> x, y;
            bool success =
                self.calculateBeta(beta_mass, beta_mol, x, y, mixCR, kij);
            return py::make_tuple(success, beta_mass, beta_mol, x, y);
          },
          py::arg("mixCR") = "CR2", py::arg("kij") = 0.0)
      .def(
          "checkSpinodal",
          [](Mixture &self, const std::string &mixCR, double kij) {
            double detQ = 0.0;
            bool inside = self.checkSpinodal(detQ, mixCR, kij);
            return py::make_tuple(inside, detQ);
          },
          py::arg("mixCR") = "CR2", py::arg("kij") = 0.0)
      .def_readwrite("speciesList", &Mixture::speciesList)
      .def_readwrite("Species_names", &Mixture::Species_names)
      .def_readwrite("Y", &Mixture::Y)
      .def_readwrite("X", &Mixture::X)
      .def_readwrite("T", &Mixture::T)
      .def_readwrite("P", &Mixture::P);
}
