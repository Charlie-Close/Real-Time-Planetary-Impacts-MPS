//
//  hdfHandler.cpp
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#include "hdfHandler.hpp"
#include <hdf5.h>
#include <iostream>

DataStruct readHDFFile(const std::string& filepath) {
    // Open the HDF5 file in read-only mode
    hid_t file_id = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error opening file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }

    // Open the "PartType0" group
    hid_t group_id = H5Gopen(file_id, "PartType0", H5P_DEFAULT);
    if (group_id < 0) {
        std::cerr << "Error opening group: GasParticles" << std::endl;
        H5Fclose(file_id);
        exit(EXIT_FAILURE);
    }

    // Helper function to read a dataset and convert to vec3 if necessary
    auto readDatasetVec3 = [&](const char* dataset_name) -> std::vector<simd_float3> {
        hid_t dataset_id = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
        if (dataset_id < 0) {
            std::cerr << "Error opening dataset: " << dataset_name << std::endl;
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Get the dataspace and dataset dimensions
        hid_t dataspace_id = H5Dget_space(dataset_id);
        hsize_t dims[2];
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
        if (dims[1] != 3) {
            std::cerr << "Dataset " << dataset_name << " does not have 3 columns as expected." << std::endl;
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Allocate space and read the data
        std::vector<float> flat_data(dims[0] * 3);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_data.data());

        // Convert to std::vector<simd_float3>
        std::vector<simd_float3> vec3_data(dims[0]);
        for (hsize_t i = 0; i < dims[0]; ++i) {
            vec3_data[i] = { flat_data[i * 3], flat_data[i * 3 + 1], flat_data[i * 3 + 2] };
        }

        // Close resources
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);

        return vec3_data;
    };

    // Helper function to read a scalar dataset (e.g., 'Pressures')
    auto readDatasetScalar = [&](const char* dataset_name) -> std::vector<float> {
        hid_t dataset_id = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
        if (dataset_id < 0) {
            std::cerr << "Error opening dataset: " << dataset_name << std::endl;
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Get the dataspace and dataset dimensions
        hid_t dataspace_id = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

        // Allocate space and read the data
        std::vector<float> data(dims[0]);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close resources
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);

        return data;
    };

    auto readDatasetInteger = [&](const char* dataset_name) -> std::vector<int> {
        hid_t dataset_id = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
        if (dataset_id < 0) {
            std::cerr << "Error opening dataset: " << dataset_name << std::endl;
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Get the dataspace and dataset dimensions
        hid_t dataspace_id = H5Dget_space(dataset_id);
        hsize_t dims[1];
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

        // Allocate space and read the data
        std::vector<int> data(dims[0]);
        H5Dread(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close resources
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);

        return data;
    };

    DataStruct data;
    data.positions = readDatasetVec3("Coordinates");
    data.densities = readDatasetScalar("Densities");
//    data.entropies = readDatasetScalar("Entropies");
    data.internalEnergy = readDatasetScalar("InternalEnergies");
    data.masses = readDatasetScalar("Masses");
    data.materialIDs = readDatasetInteger("MaterialIDs");
//    data.potentials = readDatasetScalar("Potentials");
    data.pressures = readDatasetScalar("Pressures");
    data.smoothingLengths = readDatasetScalar("SmoothingLengths");
    data.velocities = readDatasetVec3("Velocities");

    // Close the group and file
    H5Gclose(group_id);
    H5Fclose(file_id);

    return data;
}


void writeHDFFile(const std::string& filepath, const DataStruct& data) {
    // Create the HDF5 file in write mode, overwriting if it already exists
    hid_t file_id = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error creating file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }

    // Create the "PartType0" group
    hid_t group_id = H5Gcreate(file_id, "PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) {
        std::cerr << "Error creating group: PartType0" << std::endl;
        H5Fclose(file_id);
        exit(EXIT_FAILURE);
    }

    // Helper function to write a vec3 dataset
    auto writeDatasetVec3 = [&](const char* dataset_name, const std::vector<simd_float3>& vec3_data) {
        hsize_t dims[2] = { vec3_data.size(), 3 };
        hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
        hid_t dataset_id = H5Dcreate(group_id, dataset_name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        if (dataset_id < 0) {
            std::cerr << "Error creating dataset: " << dataset_name << std::endl;
            H5Sclose(dataspace_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Flatten the data
        std::vector<float> flat_data;
        flat_data.reserve(vec3_data.size() * 3);
        for (const auto& vec : vec3_data) {
            flat_data.push_back(vec.x);
            flat_data.push_back(vec.y);
            flat_data.push_back(vec.z);
        }

        // Write the data
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_data.data());

        // Close resources
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    };

    // Helper function to write a scalar dataset
    auto writeDatasetScalar = [&](const char* dataset_name, const std::vector<float>& data) {
        hsize_t dims[1] = { data.size() };
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dataset_id = H5Dcreate(group_id, dataset_name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        if (dataset_id < 0) {
            std::cerr << "Error creating dataset: " << dataset_name << std::endl;
            H5Sclose(dataspace_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Write the data
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close resources
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    };

    // Helper function to write an integer dataset
    auto writeDatasetInteger = [&](const char* dataset_name, const std::vector<int>& data) {
        hsize_t dims[1] = { data.size() };
        hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
        hid_t dataset_id = H5Dcreate(group_id, dataset_name, H5T_NATIVE_INT32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        if (dataset_id < 0) {
            std::cerr << "Error creating dataset: " << dataset_name << std::endl;
            H5Sclose(dataspace_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            exit(EXIT_FAILURE);
        }

        // Write the data
        H5Dwrite(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

        // Close resources
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
    };

    // Write datasets to the HDF5 file
    writeDatasetVec3("Coordinates", data.positions);
    writeDatasetScalar("Densities", data.densities);
    writeDatasetScalar("Entropies", data.entropies);
    writeDatasetScalar("InternalEnergies", data.internalEnergy);
    writeDatasetScalar("Masses", data.masses);
    writeDatasetInteger("MaterialIDs", data.materialIDs);
    writeDatasetScalar("Potentials", data.potentials);
    writeDatasetScalar("Pressures", data.pressures);
    writeDatasetScalar("SmoothingLengths", data.smoothingLengths);
    writeDatasetVec3("Velocities", data.velocities);
    
    std::vector<int> particleIds(data.positions.size());
    for (int i = 0; i < particleIds.size(); i++) {
        particleIds[i] = i;
    }
    
    writeDatasetInteger("ParticleIDs", particleIds);
    
    

    // Close the group and file
    H5Gclose(group_id);
    H5Fclose(file_id);
}
