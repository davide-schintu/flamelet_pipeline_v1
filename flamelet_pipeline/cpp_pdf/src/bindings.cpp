#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "integration_core.hpp"
#include "interpolation.hpp"
#include "beta_quadrature.hpp"
#include <cstring>

namespace py = pybind11;

// Wrapper function to handle Numpy arrays
py::array_t<double> integrate_py(
    py::array_t<double> Z_array,
    py::array_t<double> S_array,
    py::array_t<double> Y_array
) {
    // Request buffer info
    py::buffer_info buf_Z = Z_array.request();
    py::buffer_info buf_S = S_array.request();
    py::buffer_info buf_Y = Y_array.request();
    
    // Check dimensions
    if (buf_Z.ndim != 1) throw std::runtime_error("Z must be 1D");
    if (buf_S.ndim != 1) throw std::runtime_error("S must be 1D");
    if (buf_Y.ndim != 3) throw std::runtime_error("Y must be 3D [nP, nChi, nZ]");
    
    size_t nP = buf_Y.shape[0];
    size_t nChi = buf_Y.shape[1];
    size_t nZ = buf_Y.shape[2];
    size_t nS = buf_S.shape[0];
    
    if (size_t(buf_Z.shape[0]) != nZ) throw std::runtime_error("Z dimension mismatch");
    
    // Convert to std::vector (copying data for safety/simplicity in C++ core)
    // Note: Can optimize avoiding copy if core accepts pointers, but vector is safer for established API
    // Given zero-copy goal, copying HUGE arrays might be bad. 
    // But `integration_core` takes `vector`. 
    // For now, let's copy. 320 MB (full dataset) copy is ~0.1s. Acceptable.
    
    std::vector<double> Z_vec((double*)buf_Z.ptr, (double*)buf_Z.ptr + nZ);
    std::vector<double> S_vec((double*)buf_S.ptr, (double*)buf_S.ptr + nS);
    std::vector<double> Y_vec((double*)buf_Y.ptr, (double*)buf_Y.ptr + (nP * nChi * nZ));
    
    // Call core integration
    std::vector<double> result_vec = core::integrate_property(Z_vec, S_vec, Y_vec, nP, nChi);
    
    // Create result Numpy array
    // Shape: [nP, nChi, nZ, nS]
    std::vector<size_t> shape = {nP, nChi, nZ, nS};
    
    // Allocate numpy array
    py::array_t<double> result = py::array_t<double>(shape);
    py::buffer_info buf_out = result.request();
    
    // Copy result back
    std::memcpy(buf_out.ptr, result_vec.data(), result_vec.size() * sizeof(double));
    
    return result;
}

py::array_t<double> integrate_density_opensmoke_py(
    py::array_t<double> Z_array,
    py::array_t<double> S_array,
    py::array_t<double> rho_array
) {
    py::buffer_info buf_Z = Z_array.request();
    py::buffer_info buf_S = S_array.request();
    py::buffer_info buf_rho = rho_array.request();

    if (buf_Z.ndim != 1) throw std::runtime_error("Z must be 1D");
    if (buf_S.ndim != 1) throw std::runtime_error("S must be 1D");
    if (buf_rho.ndim != 3) throw std::runtime_error("rho must be 3D [nP, nChi, nZ]");

    size_t nP = buf_rho.shape[0];
    size_t nChi = buf_rho.shape[1];
    size_t nZ = buf_rho.shape[2];
    size_t nS = buf_S.shape[0];
    if (size_t(buf_Z.shape[0]) != nZ) throw std::runtime_error("Z dimension mismatch");

    std::vector<double> Z_vec((double*)buf_Z.ptr, (double*)buf_Z.ptr + nZ);
    std::vector<double> S_vec((double*)buf_S.ptr, (double*)buf_S.ptr + nS);
    std::vector<double> rho_vec((double*)buf_rho.ptr, (double*)buf_rho.ptr + (nP * nChi * nZ));

    std::vector<double> result_vec = core::integrate_density_opensmoke(Z_vec, S_vec, rho_vec, nP, nChi);

    py::array_t<double> result = py::array_t<double>({nP, nChi, nZ, nS});
    py::buffer_info buf_out = result.request();
    std::memcpy(buf_out.ptr, result_vec.data(), result_vec.size() * sizeof(double));
    return result;
}

py::array_t<double> integrate_favre_py(
    py::array_t<double> Z_array,
    py::array_t<double> S_array,
    py::array_t<double> rho_array,
    py::array_t<double> property_array
) {
    py::buffer_info buf_Z = Z_array.request();
    py::buffer_info buf_S = S_array.request();
    py::buffer_info buf_rho = rho_array.request();
    py::buffer_info buf_property = property_array.request();

    if (buf_Z.ndim != 1) throw std::runtime_error("Z must be 1D");
    if (buf_S.ndim != 1) throw std::runtime_error("S must be 1D");
    if (buf_rho.ndim != 3) throw std::runtime_error("rho must be 3D [nP, nChi, nZ]");
    if (buf_property.ndim != 3) throw std::runtime_error("property must be 3D [nP, nChi, nZ]");

    size_t nP = buf_rho.shape[0];
    size_t nChi = buf_rho.shape[1];
    size_t nZ = buf_rho.shape[2];
    size_t nS = buf_S.shape[0];
    if (size_t(buf_Z.shape[0]) != nZ) throw std::runtime_error("Z dimension mismatch");
    if (buf_property.shape[0] != buf_rho.shape[0] ||
        buf_property.shape[1] != buf_rho.shape[1] ||
        buf_property.shape[2] != buf_rho.shape[2]) {
        throw std::runtime_error("rho and property dimension mismatch");
    }

    std::vector<double> Z_vec((double*)buf_Z.ptr, (double*)buf_Z.ptr + nZ);
    std::vector<double> S_vec((double*)buf_S.ptr, (double*)buf_S.ptr + nS);
    std::vector<double> rho_vec((double*)buf_rho.ptr, (double*)buf_rho.ptr + (nP * nChi * nZ));
    std::vector<double> property_vec((double*)buf_property.ptr, (double*)buf_property.ptr + (nP * nChi * nZ));

    std::vector<double> rho_mean = core::integrate_density_opensmoke(Z_vec, S_vec, rho_vec, nP, nChi);
    std::vector<double> result_vec = core::integrate_favre_property(
        Z_vec, S_vec, rho_vec, property_vec, rho_mean, nP, nChi);

    py::array_t<double> result = py::array_t<double>({nP, nChi, nZ, nS});
    py::buffer_info buf_out = result.request();
    std::memcpy(buf_out.ptr, result_vec.data(), result_vec.size() * sizeof(double));
    return result;
}

py::array_t<double> integrate_1d_scalar(
    py::array_t<double> Z_array,
    py::array_t<double> Y_array,
    double Zvar
) {
    py::buffer_info buf_Z = Z_array.request();
    py::buffer_info buf_Y = Y_array.request();

    if (buf_Z.ndim != 1) throw std::runtime_error("Z must be 1D");
    if (buf_Y.ndim != 1) throw std::runtime_error("Y must be 1D");

    size_t nZ = buf_Z.shape[0];
    if (size_t(buf_Y.shape[0]) != nZ) throw std::runtime_error("Z and Y dimension mismatch");

    std::vector<double> Z_vec((double*)buf_Z.ptr, (double*)buf_Z.ptr + nZ);
    std::vector<double> Y_vec((double*)buf_Y.ptr, (double*)buf_Y.ptr + nZ);

    py::array_t<double> result = py::array_t<double>(nZ);
    py::buffer_info buf_out = result.request();
    double* out_ptr = (double*)buf_out.ptr;

    pdf::BetaQuadrature master_quad;
    pdf::FastInterpolator interp_b, interp_f, interp_center;
    interp_b.Initialize(Z_vec, master_quad.get_xb());
    interp_f.Initialize(Z_vec, master_quad.get_xf());
    interp_center.Initialize(Z_vec, master_quad.get_xcenter());

    std::vector<double> f_b(master_quad.get_xb().size());
    std::vector<double> f_f(master_quad.get_xf().size());
    std::vector<double> f_center(master_quad.get_xcenter().size());

    interp_b.Interpolate(Y_vec, f_b);
    interp_f.Interpolate(Y_vec, f_f);
    interp_center.Interpolate(Y_vec, f_center);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t iz = 0; iz < nZ; ++iz) {
        double mean = Z_vec[iz];
        double mazvard = mean * (1.0 - mean);
        double variance = Zvar;
        
        double s = (mazvard > 1e-12) ? (variance / mazvard) : 0.0;
        if (s > 1.0) s = 1.0;
        
        // Very small nonzero relative variance can introduce a visible jump when
        // switching from the identity branch to the beta-PDF branch. In practice,
        // using a first nonzero S >= 0.002 has been observed to avoid that artifact.
        if (s <= 1e-6) {
            out_ptr[iz] = Y_vec[iz];
            continue;
        }
        if (s >= 1.0 - 1e-6) {
            out_ptr[iz] = Y_vec[nZ-1]*mean + Y_vec[0]*(1.0-mean);
            continue;
        }
        if (mean < 1e-6) { out_ptr[iz] = Y_vec[0]; continue; }
        if (mean > 1.0-1e-6) { out_ptr[iz] = Y_vec[nZ-1]; continue; }

        double clamped_variance = s * mazvard;
        pdf::BetaQuadrature local_quad;
        local_quad.Set(mean, clamped_variance);
        out_ptr[iz] = local_quad.IntegralNormalized(f_b, f_f, f_center);
    }

    return result;
}

py::array_t<double> integrate_1d_array(
    py::array_t<double> Z_array,
    py::array_t<double> Y_array,
    py::array_t<double> Zvar_array
) {
    py::buffer_info buf_Z = Z_array.request();
    py::buffer_info buf_Y = Y_array.request();
    py::buffer_info buf_Zvar = Zvar_array.request();

    if (buf_Z.ndim != 1) throw std::runtime_error("Z must be 1D");
    if (buf_Y.ndim != 1) throw std::runtime_error("Y must be 1D");
    if (buf_Zvar.ndim != 1) throw std::runtime_error("Zvar must be 1D");

    size_t nZ = buf_Z.shape[0];
    if (size_t(buf_Y.shape[0]) != nZ || size_t(buf_Zvar.shape[0]) != nZ) {
        throw std::runtime_error("Dimension mismatch");
    }

    std::vector<double> Z_vec((double*)buf_Z.ptr, (double*)buf_Z.ptr + nZ);
    std::vector<double> Y_vec((double*)buf_Y.ptr, (double*)buf_Y.ptr + nZ);
    std::vector<double> Zvar_vec((double*)buf_Zvar.ptr, (double*)buf_Zvar.ptr + nZ);

    py::array_t<double> result = py::array_t<double>(nZ);
    py::buffer_info buf_out = result.request();
    double* out_ptr = (double*)buf_out.ptr;

    pdf::BetaQuadrature master_quad;
    pdf::FastInterpolator interp_b, interp_f, interp_center;
    interp_b.Initialize(Z_vec, master_quad.get_xb());
    interp_f.Initialize(Z_vec, master_quad.get_xf());
    interp_center.Initialize(Z_vec, master_quad.get_xcenter());

    std::vector<double> f_b(master_quad.get_xb().size());
    std::vector<double> f_f(master_quad.get_xf().size());
    std::vector<double> f_center(master_quad.get_xcenter().size());

    interp_b.Interpolate(Y_vec, f_b);
    interp_f.Interpolate(Y_vec, f_f);
    interp_center.Interpolate(Y_vec, f_center);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t iz = 0; iz < nZ; ++iz) {
        double mean = Z_vec[iz];
        double mazvard = mean * (1.0 - mean);
        double variance = Zvar_vec[iz];
        
        double s = (mazvard > 1e-12) ? (variance / mazvard) : 0.0;
        if (s > 1.0) s = 1.0;
        
        // Very small nonzero relative variance can introduce a visible jump when
        // switching from the identity branch to the beta-PDF branch. In practice,
        // using a first nonzero S >= 0.002 has been observed to avoid that artifact.
        if (s <= 1e-6) {
            out_ptr[iz] = Y_vec[iz];
            continue;
        }
        if (s >= 1.0 - 1e-6) {
            out_ptr[iz] = Y_vec[nZ-1]*mean + Y_vec[0]*(1.0-mean);
            continue;
        }
        if (mean < 1e-6) { out_ptr[iz] = Y_vec[0]; continue; }
        if (mean > 1.0-1e-6) { out_ptr[iz] = Y_vec[nZ-1]; continue; }

        double clamped_variance = s * mazvard;
        pdf::BetaQuadrature local_quad;
        local_quad.Set(mean, clamped_variance);
        out_ptr[iz] = local_quad.IntegralNormalized(f_b, f_f, f_center);
    }

    return result;
}

PYBIND11_MODULE(pdf_integrator_cpp, m) {
    m.doc() = "High-performance PDF Integrator using OpenSMOKE Quadrature";
    m.def("integrate", &integrate_py, "Integrate flamelet property over PDF",
          py::arg("Z"), py::arg("S"), py::arg("Y"));
    m.def("integrate_density_opensmoke", &integrate_density_opensmoke_py,
          "Integrate density with OpenSMOKE convention rho=1/<1/rho>",
          py::arg("Z"), py::arg("S"), py::arg("rho"));
    m.def("integrate_favre", &integrate_favre_py,
          "Integrate a Favre property as <rho*phi>/<rho> using OpenSMOKE density",
          py::arg("Z"), py::arg("S"), py::arg("rho"), py::arg("property"));
    m.def("integrate_1d", &integrate_1d_scalar, "Integrate 1D property over PDF with scalar Zvar",
          py::arg("Z"), py::arg("Y"), py::arg("Zvar"));
    m.def("integrate_1d", &integrate_1d_array, "Integrate 1D property over PDF with array Zvar",
          py::arg("Z"), py::arg("Y"), py::arg("Zvar"));
}
