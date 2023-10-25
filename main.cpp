#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>

#include "dpm.hpp"

int main() {
    // define the dir where everything will be saved

    std::string dir_base = "../data/tracer-sim/4/";
    tracerTest(dir_base, 4);

    return 0;
}

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