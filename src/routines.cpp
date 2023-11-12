#include <cmath>
#include <vector>
#include <algorithm>

#include "routines.hpp"
#include "dpm.hpp"
#include "disk.hpp"
#include "misc.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- Disk Functions ---------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

void adamMinimizeDiskForces(Disks2D& disks, int max_steps, double alpha, double beta1, double beta2, double epsilon) {
    // will store gradient, first moment, second moment in: force, vel, acc

    double force_norm = 0.0;
    double area = 0.0;

    // run the loop
    for (int t = 0; t < max_steps; ++t) {
        // fill the vectors with zeros
        std::fill_n(disks.force.begin(), disks.force.size(), 0.0);
        std::fill_n(disks.vel.begin(), disks.vel.size(), 0.0);
        std::fill_n(disks.acc.begin(), disks.acc.size(), 0.0);

        // calculate the forces
        disks.setInteractionForceEnergy();

        // set the force norm back to 0
        force_norm = 0.0;

        // update the moments
        for (int i = 0; i < disks.n_particles; ++i) {
            for (int dim = 0; dim < disks.simparams.N_dim; ++dim) {
                // update the first moment
                // NOTE the negative sign on the gradient is because we get the gradient of the potential from the force (which is the negative gradient of the potential)
                disks.vel[i * disks.simparams.N_dim + dim] = beta1 * disks.vel[i * disks.simparams.N_dim + dim] - (1 - beta1) * disks.force[i * disks.simparams.N_dim + dim];
                
                // update the second moment
                // NOTE - no need for the negative sign here because we are squaring the gradient
                disks.acc[i * disks.simparams.N_dim + dim] = beta2 * disks.acc[i * disks.simparams.N_dim + dim] + (1 - beta2) * disks.force[i * disks.simparams.N_dim + dim] * disks.force[i * disks.simparams.N_dim + dim];

                // calculate the bias corrected first moment
                double first_moment_bias_corrected = disks.vel[i * disks.simparams.N_dim + dim] / (1 - std::pow(beta1, t + 1));

                // calculate the bias corrected second moment
                double second_moment_bias_corrected = disks.acc[i * disks.simparams.N_dim + dim] / (1 - std::pow(beta2, t + 1));

                // update the position
                disks.pos[i * disks.simparams.N_dim + dim] -= alpha * first_moment_bias_corrected / (std::sqrt(second_moment_bias_corrected) + epsilon);

                // enforce periodic boundary conditions NOTE might not be necessary
                if (disks.pos[i * disks.simparams.N_dim + dim] > disks.simparams.box_size[dim]) {
                    disks.pos[i * disks.simparams.N_dim + dim] -= disks.simparams.box_size[dim];
                } else if (disks.pos[i * disks.simparams.N_dim + dim] < 0.0) {
                    disks.pos[i * disks.simparams.N_dim + dim] += disks.simparams.box_size[dim];
                }

                // update the force norm
                force_norm += std::pow(disks.force[i * disks.simparams.N_dim + dim], 2.0);
            }
        }

        // calculate the force norm
        force_norm = std::sqrt(force_norm);

        // check if the force norm is below the threshold
        if (force_norm < epsilon) {
            break;
        }
    }
}


// this function minimizes linear repulsive potential between disks to a desired packing fraction, iterating on the total force norm
// it is not well suited to be used for comrpessing to a nonzero pressure because the force norm is not a good indicator of the pressure
// especially bad for significant pressure
Disks2D packDisks2D_0Pressure(int num_disks, std::vector<double> diameters, std::vector<double> fraction, SimParams2D& simparams, double seed, double num_compression_steps, double compression_increment, double force_tol, double phi_target, double num_optim_steps, double learning_rate) {
    double sum = 0.0;
    int total = 0;
    for (int i = 0; i < fraction.size(); i++) {
        sum += fraction[i];
        total += 1;
    }
    if (sum != 1.0) {
        std::cout << "ERROR: fraction does not add to 1.0" << std::endl;
    }

    // make the diameter list
    std::vector<double> diams;
    for (int i = 0; i < fraction.size(); i++) {
        int num = std::round(fraction[i] * num_disks);
        for (int j = 0; j < num; j++) {
            diams.push_back(diameters[i]);
        }
    }

    // make the disks
    Disks2D disks(diams.size(), simparams);
    std::cout << "Made " << diams.size() << " disks" << std::endl;

    // randomly initialize the disk positions and diameters
    std::vector<double> rand_pos = generateRandomUniformVector(disks.n_particles * disks.simparams.N_dim, 0.0, disks.simparams.box_size[0], seed);
    for (int i = 0; i < disks.n_particles; ++i) {
        for (int dim = 0; dim < disks.simparams.N_dim; ++dim) {
            disks.pos[disks.simparams.N_dim * i + dim] = fmod(rand_pos[disks.simparams.N_dim * i + dim], disks.simparams.box_size[dim]);
        }
        disks.sigma[i] = diams[i];
    }

    // relax the disks
    for (int step = 0; step < num_compression_steps; ++step) {
        adamMinimizeDiskForces(disks, num_optim_steps, learning_rate, 0.9, 0.999, force_tol);

        // calculate the packing fraction
        // TODO need to carefully account for overlaps!
        // TODO make this more general
        double phi = 0.0;
        for (int i = 0; i < disks.n_particles; ++i) {
            phi += M_PI * std::pow(disks.sigma[i] / 2.0, 2.0);
        }
        phi /= (disks.simparams.box_size[0] * disks.simparams.box_size[1]);

        // TODO make this more general
        simparams.box_size[0] -= compression_increment;
        simparams.box_size[1] -= compression_increment;

        // check if the packing fraction is below the target
        if (phi > phi_target) {
            break;
        }

        // TODO use a generalized logger to log the data to console
    }

    return disks;
}

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- DPM Functions ----------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

void dampDpms(std::vector<DPM2D>& dpms) {
    for (int id = 0; id < dpms.size(); ++id) {
        for (int dim = 0; dim < dpms[id].simparams.N_dim; ++dim) {
            for (int i = 0; i < dpms[id].n_vertices; ++i) {
                dpms[id].force_vertex[dpms[id].simparams.N_dim * i + dim] -= dpms[id].simparams.damping_coeff * dpms[id].vel_vertex[dpms[id].simparams.N_dim * i + dim];
            }
        }
    }
}

std::vector<DPM2D> generateDpmsFromDisks(Disks2D& disks, double vertex_circumferencial_density, double radius_shrink_factor) {
    // assign the coordinates to the dpms
    std::vector<DPM2D> dpms;
    for (int i = 0; i < disks.n_particles; ++i) {
        // TODO - fix shrink factor so that we dont need the shrink radius - i think this is indicative of an issue in the dpm geometry calculation
        double radius = disks.sigma[i] / 2.0 * radius_shrink_factor;  // to be really sure there are no overlaps
        int num_vertices = std::round(vertex_circumferencial_density * 2 * M_PI * radius);
        dpms.push_back(DPM2D(disks.pos[2 * i], disks.pos[2 * i + 1], radius, num_vertices, disks.simparams, 2.0, 0.0, 0.0, 0.0, 0.0));
    }
    return dpms;
}

void shiftDpmsToVelocity(std::vector<DPM2D>& dpms, double vx, double vy) {
    std::vector<double> total_vx(dpms.size(), 0.0);
    std::vector<double> total_vy(dpms.size(), 0.0);
    double vx_mean_all = 0.0;
    double vy_mean_all = 0.0;
    int total_vertices = 0;

    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            total_vx[id] += dpms[id].vel_vertex[2 * i];
            total_vy[id] += dpms[id].vel_vertex[2 * i + 1];
        }
        total_vertices += dpms[id].n_vertices;
    }

    for (int id = 0; id < dpms.size(); ++id) {
        vx_mean_all += total_vx[id];
        vy_mean_all += total_vy[id];
    }
    vx_mean_all /= total_vertices;
    vy_mean_all /= total_vertices;

    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            dpms[id].vel_vertex[2 * i] -= vx_mean_all;
            dpms[id].vel_vertex[2 * i + 1] -= vy_mean_all;
        }
    }
}

// TODO make this more general
void zeroDpmsAngularVelocity(std::vector<DPM2D>& dpms) {
    // calculate average velocity
    for (int id = 0; id < dpms.size(); ++id) {
        double dv = 0.0;
        double com_dist_sum = 0.0;
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            double com_dist = std::sqrt(std::pow(dpms[id].pos_vertex[2 * i], 2.0) + std::pow(dpms[id].pos_vertex[2 * i + 1], 2.0));
            com_dist_sum += com_dist;
            dv += (dpms[id].vel_vertex[2 * i] * (-dpms[id].pos_vertex[2 * i + 1] / com_dist) + dpms[id].vel_vertex[2 * i + 1] * (dpms[id].pos_vertex[2 * i] / com_dist)) * com_dist;
        }
        dv /= com_dist_sum;
        // shift the dpm velocities
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            double com_dist = std::sqrt(std::pow(dpms[id].pos_vertex[2 * i], 2.0) + std::pow(dpms[id].pos_vertex[2 * i + 1], 2.0));
            dpms[id].vel_vertex[2 * i] -= dv * (-dpms[id].pos_vertex[2 * i + 1] / com_dist);
            dpms[id].vel_vertex[2 * i + 1] -= dv * (dpms[id].pos_vertex[2 * i] / com_dist);
        }
    }        
}

void scaleDpmsToTemp(std::vector<DPM2D>& dpms, double temp_target, double seed) {
    // randomly assign velocities
    int num_vertices = 0;
    for (int id = 0; id < dpms.size(); ++id) {
        num_vertices += dpms[id].n_vertices * dpms[id].simparams.N_dim;
    }
    std::vector<double> rand_vel = generateRandomNormalVector(num_vertices, 0.0, 1.0, seed);
    
    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            for (int dim = 0; dim < dpms[id].simparams.N_dim; ++dim) {
                dpms[id].vel_vertex[i * dpms[id].simparams.N_dim + dim] = rand_vel[i * dpms[id].simparams.N_dim + dim];
            }
        }
    }

    zeroDpmsAngularVelocity(dpms);
    shiftDpmsToVelocity(dpms, 0.0, 0.0);

    // calculate the temperature
    double ke = 0.0;
    num_vertices = 0;  // this can be simplified so we dont recalculate, but the num_vertices calculation would need to be changed above!
    for (int id = 0; id < dpms.size(); ++id) {
        dpms[id].calcKineticEnergies();
        ke += dpms[id].kin_eng;
        num_vertices += dpms[id].n_vertices;
    }
    double temp = 2 * ke / (dpms[0].simparams.N_dim * num_vertices * dpms[0].simparams.kb);
    // scale the velocities
    double scale_factor = std::sqrt(temp_target / temp);

    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            for (int dim = 0; dim < dpms[id].simparams.N_dim; ++dim) {
                dpms[id].vel_vertex[dpms[id].simparams.N_dim * i + dim] *= scale_factor;
            }
        }
    }
}

void verletStepDpmList(std::vector<DPM2D>& dpms) {
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].verletPositionStep();
    }

    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].innerDpmForceRoutine();
    }
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        for (int other_id = 0; other_id < dpms.size(); ++other_id) {  // TODO - test if double counting
            if (id != other_id) {
                dpms[id].setInteractionForceEnergy(dpms[other_id]);
            }
        }
    }
    if (dpms[0].simparams.damping_coeff > 0) {
        dampDpms(dpms);
    }

    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].verletVelocityStep();
    }
}

void noseHooverVelocityVerletStepDpmList(std::vector<DPM2D>& dpms, double T_target) {
    double eta_half = 0.0;
    double ke_half_sum = 0.0;
    double ke_sum = 0.0;
    
    // update the positions
    int num_particles = 0;
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].noseHooverVelocityVerletPositionStep();
        num_particles += dpms[id].n_vertices;
    }
    double K = (dpms[0].simparams.N_dim * num_particles + 1) / 2 * dpms[0].simparams.kb * T_target;

    // update the half velocities and store them in acceleration temporarily (also updates the half-eta)
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].noseHooverVelocityVerletHalfVelocityStep(ke_half_sum, ke_sum);
    }

    // update the forces
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].innerDpmForceRoutine();
    }
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        for (int other_id = 0; other_id < dpms.size(); ++other_id) {  // TODO - test if double counting
            if (id != other_id) {
                dpms[id].setInteractionForceEnergy(dpms[other_id]);
            }
        }
    }
    if (dpms[0].simparams.damping_coeff > 0) {
        dampDpms(dpms);
    }
    
    // update the eta half-step using the kinetic energy sum
    eta_half = dpms[0].simparams.eta + dpms[0].simparams.dt / (2 * dpms[0].simparams.Q) * (ke_sum - K);

    // update the eta using the half-step kinetic energy sum
    dpms[0].simparams.eta = eta_half + dpms[0].simparams.dt / (2 * dpms[0].simparams.Q) * (ke_half_sum - K);

    // update the velocities and reset the acceleration using the force
    for (int id = 0; id < dpms.size(); ++id) {  // TODO parallelize this
        dpms[id].noseHooverVelocityVerletFullVelocityStep();
    }
}

void compressDpmsNve(std::vector<DPM2D>& dpms, int max_steps, double phi_target, int compress_every, double compression_increment) {
    // TODO add the log files here

    writeMacroConsoleHeader();

    for (int step = 0; step < max_steps; ++step) {
        
        verletStepDpmList(dpms);

        // TODO log the data here

        if (step % compress_every == 0) {
            writeMacroConsoleHeader();

            // TODO make a function to calculate the packing fraction - also need to pay careful attention to DPM area calculation - seems to be underestimating
            double phi = 0.0;
            for (int id = 0; id < dpms.size(); ++id) {
                phi += dpms[id].area;
            }
            phi /= (dpms[0].simparams.box_size[0] * dpms[0].simparams.box_size[1]);

            if (phi < phi_target) {
                for (int dim = 0; dim < dpms[0].simparams.N_dim; ++dim) {
                    dpms[0].simparams.box_size[dim] -= compression_increment;
                }
            }
            else {
                break;
            }
        }
    }

    // TODO close the log files here (if needed)
}






// void tracerTest(std::string dir, int num_vertices) {
//     GeomConfig2D simparams = GeomConfig2D(10.0, 10.0, 1e-3, 1.0);
//     ForceParams simparams = ForceParams(1.0, 1.0, 1.0, 1.0, 1.0);

//     initDataFiles(dir, simparams);

//     // make the logs
//     std::ofstream macro_log = createMacroLogFile(dir + "macro_log.csv");
//     std::ofstream vertex_log = createVertexLogFile(dir + "vertex_log.csv");
//     std::ofstream dpm_log = createDpmLogFile(dir + "dpm_log.csv");

//     writeMacroConsoleHeader();

//     int num_dpms = 1;
//     int save_freq = 100;
//     int console_log_freq = 1000;
//     double T_0 = 0.0;
//     double vx_0 = 0.0;
//     double vy_0 = 0.0;
//     double cx = 5.0;
//     double cy = 5.0;
//     double R = 2.0;

//     // tracer test

//     int N_steps = 1000;

//     // make the dpm
//     DPM2D dpm = DPM2D(cx, cy, R, num_vertices, simparams, simparams, 2.0, vx_0, vy_0, T_0, 12345.0);
    
//     // make the tracer
//     double trace_radius = R + dpm.simparams.sigma / 2.0;
//     DPM2D tracer = DPM2D(cx, cy, trace_radius, 1, simparams, simparams, 2.0, vx_0, vy_0, T_0, 12345.0);
//     tracer.simparams.sigma = dpm.simparams.sigma;

//     // make a vector with two dpms
//     std::vector<DPM2D> dpms;
//     dpms.push_back(dpm);
//     dpms.push_back(tracer);

//     // loop
//     for (int i = 0; i < N_steps; ++i) {
//         // define dpms[1] (tracer) position as a function of theta
//         double theta = 2.0 * M_PI * (i) / N_steps;
//         dpms[1].pos_vertex[0] = cx + trace_radius * cos(theta);
//         dpms[1].pos_vertex[1] = cy + trace_radius * sin(theta);

//         // calculate the inner dpm forces
//         dpms[0].innerDpmForceRoutine();
//         dpms[1].force_vertex[0] = 0.0;
//         dpms[1].force_vertex[1] = 0.0;

//         // calculate the interaction forces between the dpm and the tracer
//         dpms[0].setInteractionForceEnergy(dpms[1]);

//         dpms[0].writeToLogFiles(vertex_log, dpm_log, i, 0);
//         dpms[1].writeToLogFiles(vertex_log, dpm_log, i, 1);

//         // std::cout << dpms[1].force_vertex[0] << " | " << dpms[1].force_vertex[1] << std::endl;
//     }

//     // close the logs
//     vertex_log.close();
//     dpm_log.close();
// }


// void plateCompressionSweep(std::string dir_base, int num_vertices, int N_points) {
//     // define the dir where everything will be saved
//     int iter = 0;

//     double lower_bound = std::log(1);
//     double upper_bound = std::log(500);
//     double step_size = (upper_bound - lower_bound) / (N_points - 1);

//     for (int i = 0; i < N_points; ++i) {
//         for (int j = 0; j < N_points; ++j) {
//             for (int k = 0; k < N_points; ++k) {

//                 // single dpm parallel plate compression test
//                 iter += 1;

//                 double k_a = std::floor(std::exp(lower_bound + i * step_size));
//                 double k_l = std::floor(std::exp(lower_bound + j * step_size));
//                 double k_b = std::floor(std::exp(lower_bound + k * step_size));

//                 GeomConfig2D simparams = GeomConfig2D(10.0, 10.0, 1e-3, 1.0);
//                 ForceParams simparams = ForceParams(200.0, k_l, k_b, k_a, 1.0);

//                 std::string dir = dir_base + std::to_string(iter) + "/";
//                 std::cout << dir << " of " << std::pow(N_points, 3.0) << std::endl;

//                 initDataFiles(dir, simparams);

//                 // make the logs
//                 std::ofstream macro_log = createMacroLogFile(dir + "macro_log.csv");
//                 std::ofstream vertex_log = createVertexLogFile(dir + "vertex_log.csv");
//                 std::ofstream dpm_log = createDpmLogFile(dir + "dpm_log.csv");
//                 std::ofstream config_log = createConfigLogFile(dir + "config_log.csv", simparams);

//                 writeMacroConsoleHeader();

//                 int num_steps = 20000;
//                 int save_freq = 10;
//                 int console_log_freq = 500;
//                 double R = 2.0;
//                 double dr = 0.00001;

//                 DPM2D dpm = DPM2D(5.0, 5.0, R, num_vertices, simparams, simparams, 2.0, 0.0, 0.0, 0.0, 0.0);
//                 std::vector<DPM2D> dpms;
//                 dpms.push_back(dpm);

//                 std::vector<double> wall_pos = {2.99, 7.01};
//                 // relax the dpm within the wall boundaries

//                 for (int step = 0; step < num_steps; ++step) {
//                     dpms[0].gradDescMinDpm(200000, 0.0000001, 1e-5, wall_pos);
//                     if (step % save_freq == 0) {
//                         logDpmList(dpms.size(), dpms, step, save_freq, console_log_freq, vertex_log, dpm_log, macro_log, config_log, simparams);
//                     }
//                     wall_pos[0] += dr;
//                     wall_pos[1] -= dr;
//                 }

//                 vertex_log.close();
//                 dpm_log.close();
//                 config_log.close();
//                 macro_log.close();
//             }
//         }
//     }
// }



