#ifndef DISK_HPP
#define DISK_HPP

#include <cmath>
#include <vector>

#include "sim.hpp"
#include "misc.hpp"

class Disks2D {
    public:
        std::vector<double> pos, vel, force, acc;
        std::vector<double> pos_i, pos_j, l_ij;
        std::vector<double> sigma;  // particle - particle interaction distance (width of particles)
        SimParams2D& simparams;  // sim params is inherited from the main simulation (as a reference) so that all dpms have the same parameters
        double phi, ke, pe;
        int n_particles;

        // Default Constructor
        Disks2D(int n_particles, SimParams2D& simparams);
        
        // Default Destructor
        ~Disks2D();

        // Inline Member Functions
        inline void setPos_i(const int i);
        inline void setPos_j(const int j);
        inline void verletPositionStep();
        inline void verletVelocityStep();
        inline void noseHooverVelocityVerletPositionStep();
        inline void noseHooverVelocityVerletHalfVelocityStep(double& ke_half_sum, double& ke_sum);
        inline void noseHooverVelocityVerletFullVelocityStep();
        inline void dampDisks(double damping);
        inline void langevinBAOABFullPositionStep(const double temp_target, const double gamma);
        inline void langevinBAOABFullVelocityStep(const double temp_target, const double gamma);

        // Member Functions
        void setInteractionForceEnergy();
        void setInteractionForceEnergyWCA();
        void verletStep(int step, double damping);
        void noseHooverVelocityVerletStepList(int step, double T_target, double damping);
        void shiftToVelocityMean(double vx, double vy);
        void scaleVelocitiesToTemp(double T_target, double seed);
};


// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- Inline Member Functions ------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

inline void Disks2D::setPos_i(const int i) {
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_i[dim] = this->pos[this->simparams.N_dim * i + dim];
    }
}

inline void Disks2D::setPos_j(const int j) {
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_j[dim] = this->pos[this->simparams.N_dim * j + dim];
    }
}

inline void Disks2D::verletPositionStep() {
    // update the positions
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->pos[this->simparams.N_dim * i + dim] += this->vel[this->simparams.N_dim * i + dim] * this->simparams.dt + 0.5 * this->force[this->simparams.N_dim * i + dim] * this->simparams.dt * this->simparams.dt / this->simparams.mass_vertex;
        }
    }
}

inline void Disks2D::verletVelocityStep() {
    // update the velocities and accelerations
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->vel[this->simparams.N_dim * i + dim] += 0.5 * this->simparams.dt * (this->force[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex + this->acc[this->simparams.N_dim * i + dim]);
            this->acc[this->simparams.N_dim * i + dim] = this->force[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex;
        }
    }
}

inline void Disks2D::noseHooverVelocityVerletPositionStep() {
    // update the positions
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->pos[this->simparams.N_dim * i + dim] += this->vel[this->simparams.N_dim * i + dim] * this->simparams.dt + 0.5 * (this->force[this->simparams.N_dim * i + dim] / simparams.mass_vertex - this->simparams.eta * this->vel[this->simparams.N_dim * i + dim]) * this->simparams.dt * this->simparams.dt;
        }
    }
}

inline void Disks2D::noseHooverVelocityVerletHalfVelocityStep(double& ke_half_sum, double& ke_sum) {
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            // store the half-step velocity in the acceleration for now
            this->acc[this->simparams.N_dim * i + dim] = this->vel[this->simparams.N_dim * i + dim] + 0.5 * (this->force[this->simparams.N_dim * i + dim] / simparams.mass_vertex - this->simparams.eta * this->vel[this->simparams.N_dim * i + dim]) * this->simparams.dt;
            // calculate the kinetic energy for the eta update
            ke_sum += this->vel[this->simparams.N_dim * i + dim] * this->vel[this->simparams.N_dim * i + dim] * this->simparams.mass_vertex / 2;
            // calculate the kinetic energy at the half-step for the upcoming eta update
            ke_half_sum += this->acc[this->simparams.N_dim * i + dim] * this->acc[this->simparams.N_dim * i + dim] * this->simparams.mass_vertex / 2;
        }
    }
}

inline void Disks2D::noseHooverVelocityVerletFullVelocityStep() {
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            // update the velocity using the half-step velocity (which is stored as the acceleration)
            this->vel[this->simparams.N_dim * i + dim] = (this->acc[this->simparams.N_dim * i + dim] + 0.5 * this->force[this->simparams.N_dim * i + dim] * this->simparams.dt / this->simparams.mass_vertex) / (1 + this->simparams.eta * this->simparams.dt / 2);
            // reset the acceleration using the full-step force (removes the half-step velocity from temporary storage here)
            this->acc[this->simparams.N_dim * i + dim] = this->force[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex;
        }
    }
}

inline void Disks2D::dampDisks(double damping) {
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->force[this->simparams.N_dim * i + dim] -= this->vel[this->simparams.N_dim * i + dim] * damping;
        }
    }
}

inline void Disks2D::langevinBAOABFullPositionStep(const double temp_target, const double gamma) {
    // in the process, acceleration is overwritten
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            // store the half-step velocity in the acceleration for now
            this->acc[this->simparams.N_dim * i + dim] = this->vel[this->simparams.N_dim * i + dim] + 0.5 * (this->force[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex) * this->simparams.dt;
            // overwrite the posotion with the half-step position
            this->pos_vertex[this->simparams.N_dim * i + dim] += this->simparams.dt * this->vel_vertex[this->simparams.N_dim * i + dim] * 0.5;
            // get random normal distribution with mean zero and std dev 1
            double G = generateRandomNormal();
            // update the velocity to the half-step prime velocity
            this->vel_vertex[this->simparams.N_dim * i + dim] = exp(-gamma * this->simparams.dt / this->simparams.mass_vertex) * this->acc_vertex[this->simparams.N_dim * i + dim] + std::sqrt(1 - exp(-2 * gamma * this->simparams.dt / this->simparams.mass_vertex)) * std::sqrt(this->simparams.kb * temp_target / this->simparams.mass_vertex) * G;
            // update the position to the full-step position
            this->pos_vertex[this->simparams.N_dim * i + dim] += this->simparams.dt * this->vel_vertex[this->simparams.N_dim * i + dim] * 0.5;
        }
    }
}

inline void Disks2D::langevinBAOABFullVelocityStep(const double temp_target, const double gamma) {
    for (int i = 0; i < this->n_particles; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            // update the acceleration from the force
            this->acc[this->simparams.N_dim * i + dim] = this->force[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex;
            // update the velocity to the full-step velocity
            this->vel[this->simparams.N_dim * i + dim] += 0.5 * this->simparams.dt * this->force[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex;
        }
    }
}

#endif // DISK_HPP