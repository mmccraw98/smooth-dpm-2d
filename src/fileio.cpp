#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>
#include <algorithm>

#include "H5Cpp.h"

#include "dpm.hpp"
#include "sim.hpp"
#include "fileio.hpp"

namespace fs = std::filesystem;

H5::H5File createH5File(const std::string& path) {
    try {
        // Extract the directory path from the full path
        fs::path dirPath = fs::path(path).parent_path();

        // Create directories if they do not exist
        if (!fs::exists(dirPath)) {
            if (fs::create_directories(dirPath)) {
                std::cout << "Directories created successfully.\n";
            } else {
                std::cout << "Failed to create directories.\n";
            }
        } else {
            std::cout << "Directories already exist.\n";
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << e.what() << '\n';
    }

    H5::H5File dpm_data(path, H5F_ACC_TRUNC);
    return dpm_data;
}

void writeMacroVariables(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms, int step, double pot_eng, double kin_eng, double phi, double temp, double press) {
    createAndWriteAttribute(timestepGroup, "t", step * dpms[0].simparams.dt, H5::PredType::NATIVE_INT);
    createAndWriteAttribute(timestepGroup, "pe", pot_eng, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "ke", kin_eng, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "te", pot_eng + kin_eng, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "phi", phi, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "temp", temp, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "eta", dpms[0].simparams.eta, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "press", press, H5::PredType::NATIVE_DOUBLE);
}

void writeSimParams(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms) {
    createAndWriteAttribute(timestepGroup, "num_dpms", dpms.size(), H5::PredType::NATIVE_INT);
    createAndWriteAttribute(timestepGroup, "N_dim", dpms[0].simparams.N_dim, H5::PredType::NATIVE_INT);
    std::vector<std::string> names = {"Lx", "Ly", "Lz"};
    for (int i = 0; i < dpms[0].simparams.N_dim; ++i) {
        createAndWriteAttribute(timestepGroup, names[i], dpms[0].simparams.box_size[i], H5::PredType::NATIVE_DOUBLE);
    }
    createAndWriteAttribute(timestepGroup, "kb", dpms[0].simparams.kb, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "dt", dpms[0].simparams.dt, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "k", dpms[0].simparams.k, H5::PredType::NATIVE_DOUBLE);
    createAndWriteAttribute(timestepGroup, "Q", dpms[0].simparams.Q, H5::PredType::NATIVE_DOUBLE);

    // write the potentially dpm specific parameters
    const hsize_t NUM_DPMS = dpms.size();
    const hsize_t DATA_PER_DPM = 11;  // 11 extra attributes
    hsize_t dims[2] = {NUM_DPMS, DATA_PER_DPM};   // N X (k, sigma, k_l, k_b, k_a, area_0, length_0, theta_0, mass_vertex, mass_dpm, id)
    H5::DataSpace dpm_params_dataspace(2, dims);
    H5::DataSet dpm_params_dataset = timestepGroup.createDataSet("dpm_params", H5::PredType::NATIVE_DOUBLE, dpm_params_dataspace);

    // write the data to the dataset
    double dpm_params_log[NUM_DPMS][DATA_PER_DPM];
    for (int i = 0; i < dpms.size(); ++i) {
        dpm_params_log[i][0] = dpms[i].simparams.k;
        dpm_params_log[i][1] = dpms[i].simparams.sigma;
        dpm_params_log[i][2] = dpms[i].simparams.k_l;
        dpm_params_log[i][3] = dpms[i].simparams.k_b;
        dpm_params_log[i][4] = dpms[i].simparams.k_a;
        dpm_params_log[i][5] = dpms[i].A_0;
        dpm_params_log[i][6] = dpms[i].l_0;
        dpm_params_log[i][7] = dpms[i].theta_0;
        dpm_params_log[i][8] = dpms[i].simparams.mass_vertex;
        dpm_params_log[i][9] = dpms[i].simparams.mass_dpm;
        dpm_params_log[i][10] = i;
    }
    dpm_params_dataset.write(dpm_params_log, H5::PredType::NATIVE_DOUBLE);
}

void writeDpmLevelData(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms) {
    const hsize_t NUM_DPMS = dpms.size();
    const hsize_t DATA_PER_DPM = dpms[0].simparams.N_dim * 3 + 5;  // 3 * 3 for pos, vel, force + 6 extra attributes
    hsize_t dims[2] = {NUM_DPMS, DATA_PER_DPM};   // N X ((x, y, z, vx, vy, vz, fx, fy, fz), id, area, perimeter, N_v, Pe, Ke)
    H5::DataSpace dpm_dataspace(2, dims);
    H5::DataSet dpm_dataset = timestepGroup.createDataSet("dpm_data", H5::PredType::NATIVE_DOUBLE, dpm_dataspace);

    // write the data to the dataset
    double dpm_data_log[NUM_DPMS][DATA_PER_DPM];
    for (int i = 0; i < dpms.size(); ++i) {
        for (int dim = 0; dim < dpms[i].simparams.N_dim; ++dim) {
            dpm_data_log[i][dim] = dpms[i].pos_dpm[dim];
            dpm_data_log[i][dim + dpms[i].simparams.N_dim] = dpms[i].vel_dpm[dim];
            dpm_data_log[i][dim + 2 * dpms[i].simparams.N_dim] = dpms[i].force_dpm[dim];
        }
        dpm_data_log[i][3 * dpms[i].simparams.N_dim] = i;
        dpm_data_log[i][3 * dpms[i].simparams.N_dim + 1] = dpms[i].area;
        dpm_data_log[i][3 * dpms[i].simparams.N_dim + 2] = dpms[i].perimeter;
        dpm_data_log[i][3 * dpms[i].simparams.N_dim + 3] = dpms[i].n_vertices;
        dpm_data_log[i][3 * dpms[i].simparams.N_dim + 4] = dpms[i].pot_eng;
        dpm_data_log[i][3 * dpms[i].simparams.N_dim + 5] = dpms[i].kin_eng;
    }
    dpm_dataset.write(dpm_data_log, H5::PredType::NATIVE_DOUBLE);
}

void writeVertexLevelData(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms, int num_vertices) {
    const hsize_t NUM_VERTICES = num_vertices;
    const hsize_t DATA_PER_VERTEX = dpms[0].simparams.N_dim * 3 + 2;  // 3 * 3 for pos, vel, force + 2 extra attributes
    hsize_t dims[2] = {NUM_VERTICES, DATA_PER_VERTEX};   // N X ((x, y, z, vx, vy, vz, fx, fy, fz), id, dpm_id)
    H5::DataSpace vertex_dataspace(2, dims);
    H5::DataSet vertex_dataset = timestepGroup.createDataSet("vertex_data", H5::PredType::NATIVE_DOUBLE, vertex_dataspace);

    // write the data to the dataset
    double vertex_data_log[num_vertices][DATA_PER_VERTEX];
    int vertex_index = 0;
    for (int i = 0; i < dpms.size(); ++i) {
        for (int j = 0; j < dpms[i].n_vertices; ++j) {
            for (int dim = 0; dim < dpms[i].simparams.N_dim; ++dim) {
                vertex_data_log[vertex_index][dim] = dpms[i].pos_vertex[j * dpms[i].simparams.N_dim + dim];
                vertex_data_log[vertex_index][dim + dpms[i].simparams.N_dim] = dpms[i].vel_vertex[j * dpms[i].simparams.N_dim + dim];
                vertex_data_log[vertex_index][dim + 2 * dpms[i].simparams.N_dim] = dpms[i].force_vertex[j * dpms[i].simparams.N_dim + dim];
            }
            vertex_data_log[vertex_index][3 * dpms[i].simparams.N_dim] = j;
            vertex_data_log[vertex_index][3 * dpms[i].simparams.N_dim + 1] = i;
            vertex_index += 1;
        }
    }
    vertex_dataset.write(vertex_data_log, H5::PredType::NATIVE_DOUBLE);
}

void writeData(H5::H5File& dpm_data, const std::vector<DPM2D>& dpms, bool save_vertex, bool save_params, int step, double pot_eng, double kin_eng, double phi, double temp, double press, int num_vertices) {

    // make a group for the current time step
    std::string groupName = "step_" + std::to_string(step);
    H5::Group timestepGroup = dpm_data.createGroup(groupName);
    writeMacroVariables(timestepGroup, dpms, step, pot_eng, kin_eng, phi, temp, press);

    // write the simulation parameters
    if (save_params) {
        writeSimParams(timestepGroup, dpms);
    }

    // write the dpm data
    writeDpmLevelData(timestepGroup, dpms);

    // write the vertex data
    if (save_vertex) {
        writeVertexLevelData(timestepGroup, dpms, num_vertices);
    }
}


void writeMacroConsoleHeader() {  // TODO make this more general
    std::cout << std::string(12 * 8 + 6, '_') << std::endl;
    std::cout << std::setw(12) << "step" << " | " 
              << std::setw(12) << "time" << " | " 
              << std::setw(12) << "pe" << " | " 
              << std::setw(12) << "ke" << " | " 
              << std::setw(12) << "te" << " | "
              << std::setw(12) << "phi" << " | " 
              << std::setw(12) << "temp" << "\n";
    std::cout << std::string(12 * 8 + 6, '_') << std::endl;
}

void writeMacroConsoleLine(int step, std::vector<DPM2D>& dpms, double pe, double ke, double phi, double temp) {
    std::cout << std::setw(12) << std::fixed << std::setprecision(3) << step << " | " 
    << std::setw(12) << std::scientific << std::setprecision(3) << (dpms[0].simparams.dt * step) << " | " 
    << std::setw(12) << std::scientific << std::setprecision(3) << pe << " | " 
    << std::setw(12) << std::scientific << std::setprecision(3) << ke << " | " 
    << std::setw(12) << std::scientific << std::setprecision(3) << ke + pe << " | " 
    << std::setw(12) << std::scientific << std::setprecision(3) << phi << " | " 
    << std::setw(12) << std::scientific << std::setprecision(3) << temp << "\n";
}

void logDpms(std::vector<DPM2D>& dpms, H5::H5File& dpm_data, int step, int log_every, int console_log_every, int rewrite_header_every, bool save_vertex, bool save_params) {
    // calculate the macroscopic variables
    double pot_eng = 0.0;
    double kin_eng = 0.0;
    double phi = 0.0;
    int num_vertices = 0;
    double press = 0.0;  // TODO calculate pressure
    for (int id = 0; id < dpms.size(); ++id) {
        pot_eng += dpms[id].pot_eng;
        kin_eng += dpms[id].kin_eng;
        phi += dpms[id].area;
        num_vertices += dpms[id].n_vertices;
    }
    double temp = 2 * kin_eng / (dpms[0].simparams.N_dim * num_vertices * dpms[0].simparams.kb);

    if (step % rewrite_header_every == 0) {
        writeMacroConsoleHeader();
    }

    if (step % console_log_every == 0) {
        writeMacroConsoleLine(step, dpms, pot_eng, kin_eng, phi, temp);
    }

    if (step % log_every == 0) {
        writeData(dpm_data, dpms, save_vertex, save_params, step, pot_eng, kin_eng, phi, temp, press, num_vertices);
    }
}