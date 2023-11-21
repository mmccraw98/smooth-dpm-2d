#include <vector>
#include <cmath>

#include "disk.hpp"
#include "sim.hpp"
#include "misc.hpp"

Disks2D::Disks2D(int n_particles, SimParams2D& simparams) : simparams(simparams) {
    // Default Constructor
    this->pos = std::vector<double>(n_particles * simparams.N_dim, 0.0);
    this->pos_i = std::vector<double>(simparams.N_dim, 0.0);
    this->pos_j = std::vector<double>(simparams.N_dim, 0.0);
    this->l_ij = std::vector<double>(simparams.N_dim, 0.0);
    this->vel = std::vector<double>(n_particles * simparams.N_dim, 0.0);
    this->force = std::vector<double>(n_particles * simparams.N_dim, 0.0);
    this->acc = std::vector<double>(n_particles * simparams.N_dim, 0.0);
    this->sigma = std::vector<double>(n_particles, 0.0);
    this->simparams = simparams;
    this->phi = 0.0;
    this->ke = 0.0;
    this->pe = 0.0;
    this->n_particles = n_particles;
}

Disks2D::~Disks2D() {
    // Default Destructor
}

// TODO - unify these functions to work with a defined potential
void Disks2D::setInteractionForceEnergy() {
    // loop over all particles
    for (int i = 0; i < this->n_particles; ++i) {
        this->setPos_i(i);
        // loop over other particles
        for (int j = i + 1; j < this->n_particles; ++j) {
            this->setPos_j(j);
            setDistVect(this->l_ij, this->pos_i, this->pos_j, this->simparams);
            double dist = getVectLength(this->l_ij, this->simparams);
            double magnitude = (dist - (this->sigma[i] + this->sigma[j]) / 2.0);
            if (dist < (this->sigma[i] + this->sigma[j]) / 2.0) {
                for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
                    this->force[this->simparams.N_dim * i + dim] -= this->simparams.k * magnitude * this->l_ij[dim] / dist;
                    this->force[this->simparams.N_dim * j + dim] += this->simparams.k * magnitude * this->l_ij[dim] / dist;
                }
                this->pe += this->simparams.k * magnitude * magnitude;
            }
        }
    }
}

void Disks2D::setInteractionForceEnergyWCA() {
    // loop over all particles
    for (int i = 0; i < this->n_particles; ++i) {
        this->setPos_i(i);
        // loop over other particles
        for (int j = i + 1; j < this->n_particles; ++j) {
            this->setPos_j(j);
            setDistVect(this->l_ij, this->pos_i, this->pos_j, this->simparams);
            double dist = getVectLength(this->l_ij, this->simparams);
            double sigma_eff = (this->sigma[i] + this->sigma[j]) / 2.0;
            if (dist < sigma_eff * 1.122462048309373) {

                double sigma12 = std::pow(sigma_eff, 12) / std::pow(dist, 12);
                double sigma6 = std::pow(sigma_eff, 6) / std::pow(dist, 6);
                double force_magnitude = 24 * sigma_eff * (2 * sigma12 - sigma6) / dist;
                double energy_magnitude = 4 * this->simparams.k * (sigma12 - sigma6) + this->simparams.k;

                for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
                    this->force[this->simparams.N_dim * i + dim] += force_magnitude * this->l_ij[dim] / dist;
                    this->force[this->simparams.N_dim * j + dim] -= force_magnitude * this->l_ij[dim] / dist;
                }
                this->pe += energy_magnitude;
            }
        }
    }
}

void Disks2D::verletStep(int step, double damping) {
    this->verletPositionStep();
    this->setInteractionForceEnergy();
    if (damping > 0) {
        this->dampDisks(damping);
    }
    this->verletVelocityStep();
}

void Disks2D::noseHooverVelocityVerletStepList(int step, double T_target, double damping) {
    double eta_half = 0.0;
    double ke_half_sum = 0.0;
    double ke_sum = 0.0;
    
    // update the positions
    this->verletPositionStep();
    this->noseHooverVelocityVerletPositionStep();
    double K = (this->simparams.N_dim * this->n_particles + 1) / 2 * this->simparams.kb * T_target;

    // update the half velocities and store them in acceleration temporarily (also updates the half-eta)
    this->noseHooverVelocityVerletHalfVelocityStep(ke_half_sum, ke_sum);

    // update the forces
    this->setInteractionForceEnergy();
    if (damping > 0) {
        this->dampDisks(damping);
    }
    
    // update the eta half-step using the kinetic energy sum
    eta_half = this->simparams.eta + this->simparams.dt / (2 * this->simparams.Q) * (ke_sum - K);

    // update the eta using the half-step kinetic energy sum
    this->simparams.eta = eta_half + this->simparams.dt / (2 * this->simparams.Q) * (ke_half_sum - K);

    // update the velocities and reset the acceleration using the force
    this->noseHooverVelocityVerletFullVelocityStep();
}

void Disks2D::shiftToVelocityMean(double vx, double vy) {
    double vx_mean = 0.0;
    double vy_mean = 0.0;
    for (int i = 0; i < this->n_particles; ++i) {
        vx_mean += this->vel[this->simparams.N_dim * i];
        vy_mean += this->vel[this->simparams.N_dim * i + 1];
    }
    vx_mean /= this->n_particles;
    vy_mean /= this->n_particles;

    for (int i = 0; i < this->n_particles; ++i) {
        this->vel[this->simparams.N_dim * i] -= vx_mean - vx;
        this->vel[this->simparams.N_dim * i + 1] -= vy_mean - vy;
    }
}

void Disks2D::scaleVelocitiesToTemp(double T_target, double seed) {
    // randomly assign velocities

    std::vector<double> rand_vel = generateRandomNormalVector(this->n_particles * this->simparams.N_dim, 0.0, 1.0, seed);
    
    // assign the velocities
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->vel[this->simparams.N_dim * i + dim] = rand_vel[this->simparams.N_dim * i + dim];
        }
    }

    // zero the center of mass velocity
    this->shiftToVelocityMean(0.0, 0.0);

    // calculate the temperature
    double ke = 0.0;
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            ke += this->vel[this->simparams.N_dim * i + dim] * this->vel[this->simparams.N_dim * i + dim] * this->simparams.mass_vertex / 2;
        }
    }
    double temp = 2 * ke / (this->simparams.N_dim * this->n_particles * this->simparams.kb);
    // scale the velocities
    double scale_factor = std::sqrt(T_target / temp);

    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->vel[this->simparams.N_dim * i + dim] *= scale_factor;
        }
    }
}
