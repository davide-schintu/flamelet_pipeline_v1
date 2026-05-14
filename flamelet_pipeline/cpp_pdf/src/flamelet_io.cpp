#include "flamelet_io.hpp"
#include <iostream>
#include <stdexcept>

namespace io {

FlameletData readInputH5(
    const std::string& filename,
    const std::vector<std::string>& property_names
) {
    FlameletData data;
    data.property_names = property_names;
    
    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        
        // Read Z grid from root
        {
            H5::DataSet dset = file.openDataSet("Z");
            H5::DataSpace space = dset.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            data.nZ = dims[0];
            data.Z_grid.resize(data.nZ);
            dset.read(data.Z_grid.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Read pressure grid
        if (H5Lexists(file.getId(), "pressure", H5P_DEFAULT) > 0) {
            H5::DataSet dset = file.openDataSet("pressure");
            H5::DataSpace space = dset.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            data.nP = dims[0];
            data.pressure.resize(data.nP);
            dset.read(data.pressure.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Read chi grid
        if (H5Lexists(file.getId(), "chi", H5P_DEFAULT) > 0) {
            H5::DataSet dset = file.openDataSet("chi");
            H5::DataSpace space = dset.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            data.nChi = dims[0];
            data.chi.resize(data.nChi);
            dset.read(data.chi.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Open data group
        H5::Group dataGroup = file.openGroup("data");
        
        // Read a reference dataset to get dimensions (use T)
        {
            H5::DataSet dset = dataGroup.openDataSet("T");
            H5::DataSpace space = dset.getSpace();
            int ndims = space.getSimpleExtentNdims();
            
            if (ndims == 3) {
                hsize_t dims[3];
                space.getSimpleExtentDims(dims);
                data.nP = dims[0];
                data.nChi = dims[1];
                data.nZ = dims[2];
            } else {
                throw std::runtime_error("Expected 3D dataset for T");
            }
        }
        
        // Read each property
        data.input_data.resize(property_names.size());
        size_t total_size = data.nP * data.nChi * data.nZ;
        
        for (size_t i = 0; i < property_names.size(); ++i) {
            const auto& prop = property_names[i];
            
            if (H5Lexists(dataGroup.getId(), prop.c_str(), H5P_DEFAULT) <= 0) {
                std::cerr << "Warning: Property '" << prop << "' not found, skipping.\n";
                data.input_data[i].resize(total_size, 0.0);
                continue;
            }
            
            H5::DataSet dset = dataGroup.openDataSet(prop);
            data.input_data[i].resize(total_size);
            dset.read(data.input_data[i].data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        file.close();
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("HDF5 error reading " + filename + ": " + e.getDetailMsg());
    }
    
    return data;
}

void writeOutputH5(
    const std::string& filename,
    const FlameletData& data
) {
    try {
        // Open existing file for read/write
        H5::H5File file(filename, H5F_ACC_RDWR);
        
        // Open or create data group
        H5::Group dataGroup;
        if (H5Lexists(file.getId(), "data", H5P_DEFAULT) > 0) {
            dataGroup = file.openGroup("data");
        } else {
            dataGroup = file.createGroup("data");
        }
        
        // Write each property
        hsize_t dims4[4] = {data.nP, data.nChi, data.nZ, data.nS};
        H5::DataSpace space4(4, dims4);
        
        for (size_t i = 0; i < data.property_names.size(); ++i) {
            const auto& prop = data.property_names[i];
            
            // Delete existing dataset if present
            if (H5Lexists(dataGroup.getId(), prop.c_str(), H5P_DEFAULT) > 0) {
                H5Ldelete(dataGroup.getId(), prop.c_str(), H5P_DEFAULT);
            }
            
            H5::DataSet dset = dataGroup.createDataSet(
                prop, H5::PredType::NATIVE_DOUBLE, space4
            );
            dset.write(data.output_data[i].data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        file.close();
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("HDF5 error writing " + filename + ": " + e.getDetailMsg());
    }
}

void copyBaseStructure(
    const std::string& input_file,
    const std::string& output_file,
    const std::vector<double>& S_grid
) {
    try {
        H5::H5File fin(input_file, H5F_ACC_RDONLY);
        H5::H5File fout(output_file, H5F_ACC_TRUNC);  // Create new
        
        // Copy Z grid
        {
            H5::DataSet dset_in = fin.openDataSet("Z");
            H5::DataSpace space = dset_in.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            
            std::vector<double> Z(dims[0]);
            dset_in.read(Z.data(), H5::PredType::NATIVE_DOUBLE);
            
            H5::DataSet dset_out = fout.createDataSet("Z", H5::PredType::NATIVE_DOUBLE, space);
            dset_out.write(Z.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Copy pressure if exists
        if (H5Lexists(fin.getId(), "pressure", H5P_DEFAULT) > 0) {
            H5::DataSet dset_in = fin.openDataSet("pressure");
            H5::DataSpace space = dset_in.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            
            std::vector<double> p(dims[0]);
            dset_in.read(p.data(), H5::PredType::NATIVE_DOUBLE);
            
            H5::DataSet dset_out = fout.createDataSet("pressure", H5::PredType::NATIVE_DOUBLE, space);
            dset_out.write(p.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Copy chi if exists
        if (H5Lexists(fin.getId(), "chi", H5P_DEFAULT) > 0) {
            H5::DataSet dset_in = fin.openDataSet("chi");
            H5::DataSpace space = dset_in.getSpace();
            hsize_t dims[1];
            space.getSimpleExtentDims(dims);
            
            std::vector<double> chi(dims[0]);
            dset_in.read(chi.data(), H5::PredType::NATIVE_DOUBLE);
            
            H5::DataSet dset_out = fout.createDataSet("chi", H5::PredType::NATIVE_DOUBLE, space);
            dset_out.write(chi.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Write S grid (Zvar_norm)
        {
            hsize_t dims[1] = {S_grid.size()};
            H5::DataSpace space(1, dims);
            H5::DataSet dset = fout.createDataSet("Zvar_norm", H5::PredType::NATIVE_DOUBLE, space);
            dset.write(S_grid.data(), H5::PredType::NATIVE_DOUBLE);
        }
        
        // Create data group
        fout.createGroup("data");
        
        fin.close();
        fout.close();
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("HDF5 error in copyBaseStructure: " + e.getDetailMsg());
    }
}

} // namespace io
