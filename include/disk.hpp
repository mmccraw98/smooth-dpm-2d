#ifndef DISK_HPP
#define DISK_HPP

#include <cmath>
#include <vector>

#include "sim.hpp"

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


#endif // DISK_HPP