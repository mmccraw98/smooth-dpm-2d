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
#include "disk.hpp"
#include "routines.hpp"
#include "fileio.hpp"


// TODO: hdf5 loading
// TODO: stress / strain / bulk modulus, etc.
// TODO: rescale units


int main() {
    double seed = 42;
    double temp_target = 1.0;
    int num_dpms = 10;
    double phi_target = 0.7;
    int num_steps = 10000;

    int log_every = 100;
    int console_log_every = 1000;
    int rewrite_header_every = 10000;
    bool save_vertex = true;
    bool save_params = true;

    // define the simulation parameters
    SimParams2D simparams = SimParams2D(20.0, 20.0, 1e-3, 1.0, 1.0, 100.0, 100.0, 100.0, 1.0);

    // make the disk packing
    Disks2D disks = packDisks2D_0Pressure(num_dpms, {1.0, 1.4}, {0.5, 0.5}, simparams, seed, 1000, 0.1, 1e-5, phi_target, 1000, 0.1);

    // make the dpms
    std::vector<DPM2D> dpms = generateDpmsFromDisks(disks, 2.0, 0.8);
    scaleDpmsToTemp(dpms, temp_target, seed);

    // create the file
    H5::H5File dpm_data = createH5File("./data/test.h5");

    // KINETIC ENERGY IS NOT BEING UPDATED CORRECTLY

    // run the dynamics
    for (int step = 0; step < num_steps; ++step) {
        noseHooverVelocityVerletStepDpmList(dpms, temp_target);

        if ((step % log_every == 0) || (step % console_log_every == 0) || (step % rewrite_header_every == 0)) {
            logDpms(dpms, dpm_data, step, log_every, console_log_every, rewrite_header_every, save_vertex, save_params);
            if (step > 0) {
                save_params = false;
            }
        }
    }

    // once done with simulation, close the file
    dpm_data.close();

    // int num_steps = 1000;
    // int save_freq = 100;
    // int console_log_freq = 1000;
    // double eta = 1.0;
    // double Q = 1.0;
    // double temp_target = 5.0;
    // double phi_target = 0.8;
    // int num_dpms = 10;

    // // define the directory where everything will be saved
    // std::string dir = "./data/test/";

    // GeomConfig2D geomconfig = GeomConfig2D(20.0, 20.0, 1e-4, 1.0);
    // ForceParams forceparams = ForceParams(1.0, 100.0, 100.0, 100.0, 1.0);
    // forceparams.eta = eta;
    // forceparams.Q = Q;

    // // make the dpms using a disk packing
    // std::vector<double> radii = {1.0, 1.4};
    // std::vector<double> fraction = {0.5, 0.5};
    // double tol = 1e-5;
    // double seed = 42;
    // PosRadius pos_rad = generateDiskPackCoords(num_dpms, radii, fraction, geomconfig, forceparams, seed, 0.1, tol, phi_target, 1000);
    // std::vector<DPM2D> dpms = generateDpmsFromDiskPack(pos_rad, geomconfig, forceparams, 3, std::pow(2, -1 / 6) * 0.87);

    // // scale the dpm velocities to be at a specified temperature
    // scaleDpmsToTemp(dpms, geomconfig, forceparams, temp_target, seed);

    // // define the dir where everything will be saved
    // initDataFiles(dir, geomconfig);

    // // make the logs
    // std::ofstream macro_log = createMacroLogFile(dir + "macro_log.csv");
    // std::ofstream vertex_log = createVertexLogFile(dir + "vertex_log.csv");
    // std::ofstream dpm_log = createDpmLogFile(dir + "dpm_log.csv");
    // std::ofstream config_log = createConfigLogFile(dir + "config_log.csv", geomconfig);

    // // assign the coordinates to the dpms

    // // initialize the quasistatic packing simulation in the normal way
    // writeMacroConsoleHeader();

    // for (int step = 0; step < num_steps; ++step) {
    //     // verletStepDpmList(dpms.size(), dpms, step, 0.0);
    //     noseHooverVelocityVerletStepDpmList(dpms.size(), dpms, step, eta, temp_target, Q, 0.0);

    //     // write the data to the files
    //     if (step % save_freq == 0) {
    //         dpms[0].forceparams.eta = eta;
    //         logDpmList(dpms.size(), dpms, step, save_freq, console_log_freq, vertex_log, dpm_log, macro_log, config_log, geomconfig);
    //     }
    // }

    // vertex_log.close();
    // dpm_log.close();
    // config_log.close();
    // macro_log.close();


    return 0;
}

// fix loading script
// make simulation script stop if nan detected when saving
// have option to restart simulation from prior step if nan detected (try to restart up to N times)

// bulk modulus and stuff

// dpm area calculation may be very innacurate - should use a better method - phi seems to be far undercalculated

// should calculate neighboring dpm non-bonded repulsion - should not get into scenarios where dpm self intersects with opposite-ended segments

// Note: if the timestep is too large, or the energy is too high, or the perimeter stiffness is too small, the dpm vertices can pass through each other
// which gives a very large increase in bending energy (the angle changes)

// figure out the units

// make logs separate files every N log steps

// split the code into different files: dpm, geomconfig, forceparams, fileio, etc.
// make a class for the simulation
// make geomconfig more distinct
// denote which functions only work in 2D

// should probably unify the compression and simulation codes

// pressure calculation

// pack, make energy conservation plot

// code for histogram

// add dpm-level variable tracking the number of contacts

// account for wca sigma better in dpm size calculation

// clean up the wca force

// better dpm size calculation

// TODO:
// zero out velocity for list of dpms
// zero out angular velocity for list of dpms
// set velocity for list of dpms given temperature
// fix the nvt
// load restart file

// code clean up:
// consolidate code (pass an energy function)
// clean up energy calculation - switch it to just PE

// the segments are effectively another type of particle - so the scaling of computation is bad since N->2N with segments

// there is no reason to pack dpms to a zero pressure state - in that case, just pack using disks
// dpms only offer an advantage when there is internal pressure leading to deformation

// better area calculation for dpm (including overlap)

// remove all dependencies on num_particles when it could be easily calculated from dpms.size()

// maybe should be compressing to a desired pressure? idk ask arthur again

// add andrews correction factor (could maybe use thte theta_ij value to determine if concave)
// validate energy conservation with and without correction factor
// nvt simulation
// fix compression algorithm so that it just packs to either desired phi or jamming point - with bidispersity
// compression with dpms
// orientation field, velocity field, msd, isf, rdf, voronoi, energy, pressure, bulk modulus, etc.
// publish python package: visualisation, analysis, etc. (make classes for things and clean up code)
// clean up c++ code, make comments, make more classes

// generalize the code to call a geomconfig object and subclass it to get a geomconfig2d object

// slides:
// disk tracer and reduction to disk-disk interaction
// compression valdiation, how to obtain cylinder contact from dpm - READ PAPERS IF ANY
// movie of packed dpm

// dpm packing - bidisperse!

// calculate bulk modulus - plot bulk modulus vs pressure - the slope gives stiffness

// add WCA interaction (repulsive lennard jones)

// begin glass formation runs

// analysis code: msd, isf, rdf, voronoi, energy, pressure, bulk modulus, etc.

// make more efficient animator (i.e. make all plot objects and just update the plot object data)
// compression of dpm between two walls
// hertz contact needs to be quasistatic i.e. will be done with energy minimization
// pack with bidispersity
// linear repulsion seems to allow pass through
// test interaction force loop double counting
// initialize velocities of multi dpm system
// github
// multiple dpms
// comments
// energy minimization
// nvt modification
// read restart file
// the area should include sigma / 2
// the packing fraction should remove overlaps
// neighbor list
// parallelization
// NOTE can parallelize the particle - dpm interaction but only between particle i and all vertices in dpm, 
// thrust

// TODO:
// automatically generate sets of dpms
// instead of taking the values, making stack arrays, and setting them back in, just directly modify the values by index
// track rotational angle of the dpm
// track moment of inertia of the dpm
// energy conservation test
// energy minimisation script
// build in performance mode
// do we use std::fill_n or just a for loop to zero arrays? - test performance
// parallelise the code
// think about memory pooling in the future
// move to thrust for gpu acceleration

// ----------- NOTES FROM PYTHON SCRIPTS:

// # is dpm equivalent to nematic liquid crystal? (i.e. is there a director field?)
// # very interesting result here:
// # if we plot the average major and minor axis lengths as functions of time, we see that there is a near immediate and permanent transition from even axis lengths to
// # large differences in the axes.  this shows that the dpms seem to prefer to be ellipsoidal rather than circular.  here is a sample code:
// # dpm_df, vertex_df, config_df, macro_df = load_sim_dataframes('sim1')
// # append_major_minor_axes_to_dpmdf(dpm_df, vertex_df)
// # plt.plot(dpm_df.groupby('step')['major_axis'].mean())
// # plt.plot(dpm_df.groupby('step')['minor_axis'].mean())
// # ######
// # this can be potentially understood if ellipses pack at a higher packing fraction than circles


// # calcualte velocity field

// # calculate msd (maybe need it to be more general using lag times?)
// # calculate general time correlator of a quantity
// # calculate general space correlator of a quantity
// # calculate general space-time correlator of a quantity

// coordination number

// # calculate rdf

// # calculate isf

// # calculate shape parameter

// # number of neighbors

// # dynamic matrix

// # contact network

// # static structure factor

// # calculate order parameter

// # x4 dynamic susceptibility

// # collective intermediate scattering function

// # angell plot

// # (these all have scipy functions for calculation AND plotting!)
// # calcualte voronoi
// # calculate delaunay
// # calculate convex hull