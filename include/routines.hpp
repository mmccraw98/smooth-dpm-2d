#ifndef ROUTINES_HPP
#define ROUTINES_HPP

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

#include "disk.hpp"
#include "sim.hpp"
#include "dpm.hpp"
#include "misc.hpp"
#include "fileio.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- Disk Functions ---------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

void adamMinimizeDiskForces(Disks2D& disks, int max_steps, double alpha, double beta1, double beta2, double epsilon);
Disks2D packDisks2D_0Pressure(int num_disks, std::vector<double> diameters, std::vector<double> fraction, SimParams2D& simparams, double seed, double num_compression_steps, double compression_increment, double force_tol, double phi_target, double num_optim_steps, double learning_rate);

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- DPM Functions ----------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

void dampDpms(std::vector<DPM2D>& dpms);
std::vector<DPM2D> generateDpmsFromDisks(Disks2D& disks, double vertex_circumferencial_density, double radius_shrink_factor);
void shiftDpmsToVelocity(std::vector<DPM2D>& dpms, double vx, double vy);
void zeroDpmsAngularVelocity(std::vector<DPM2D>& dpms);
void scaleDpmsToTemp(std::vector<DPM2D>& dpms, double temp_target, double seed);
void verletStepDpmList(std::vector<DPM2D>& dpms);
void noseHooverVelocityVerletStepDpmList(std::vector<DPM2D>& dpms, double T_target);
void compressDpmsNve(std::vector<DPM2D>& dpms, int max_steps, double phi_target, int compress_every, double compression_increment);

#endif // ROUTINES_HPP
