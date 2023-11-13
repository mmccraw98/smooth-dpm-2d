#ifndef FILEIO_HPP
#define FILEIO_HPP

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

namespace fs = std::filesystem;

template <typename T>
void createAndWriteAttribute(H5::Group& group, const std::string& attributeName, const T& value, const H5::PredType& predType) {
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName, predType, attr_dataspace);
    attribute.write(predType, &value);
}

template <typename T>
void readAttribute(H5::H5Object& obj, const std::string& attributeName, T& value) {
    H5::Attribute attribute = obj.openAttribute(attributeName);
    H5::DataType dataType = attribute.getDataType();
    attribute.read(dataType, &value);
}

H5::H5File createH5File(const std::string& path);
void writeMacroVariables(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms, int step, double pot_eng, double kin_eng, double phi, double temp, double press);
void writeSimParams(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms);
void writeDpmLevelData(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms);
void writeVertexLevelData(H5::Group& timestepGroup, const std::vector<DPM2D>& dpms, int num_vertices);
void writeData(H5::H5File& dpm_data, const std::vector<DPM2D>& dpms, bool save_vertex, bool save_params, int step, double pot_eng, double kin_eng, double phi, double temp, double press, int num_vertices);
void writeMacroConsoleHeader();
void writeMacroConsoleLine(int step, std::vector<DPM2D>& dpms, double pe, double ke, double phi, double temp);
void logDpms(std::vector<DPM2D>& dpms, H5::H5File& dpm_data, int step, int log_every, int console_log_every, int rewrite_header_every, bool save_vertex, bool save_params);
void readSimParams(H5::H5File& dpm_data, int& latest_sim_params_update_step, SimParams2D& simparams, int& num_dpms);
void readDpmData(H5::H5File& dpm_data, const int desired_step, std::vector<DPM2D>& dpms, const int latest_sim_params_update_step, SimParams2D& simparams);
void readDpmDataFromStep(std::string dpm_data_path, int &desired_step, std::vector<DPM2D>& dpms, SimParams2D& simparams);
#endif // FILEIO_HPP