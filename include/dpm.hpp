#ifndef DPM_HPP
#define DPM_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <sstream>
#include <filesystem>

#include "sim.hpp"
#include "misc.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- DPM Class --------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

class DPM2D {
    public:
        // Member variables
        std::vector<double> pos_vertex, vel_vertex, force_vertex, acc_vertex;
        std::vector<double> bond_lengths, bond_angles;
        std::vector<double> pos_dpm, vel_dpm, force_dpm;
        std::vector<double> pos_i, pos_j, pos_k;
        std::vector<double> l_ij, l_ki;
        std::vector<bool> vertex_is_active;
        double A_0, l_0, theta_0;
        double pot_eng, kin_eng;
        int n_vertices;
        double sigma;
        double area, perimeter;
        SimParams2D& simparams;  // sim params is inherited from the main simulation (as a reference) so that all dpms have the same parameters

        // Inline Member Functions
        void initializeVectors();
        inline void verletPositionStep();
        inline void verletVelocityStep();
        inline void noseHooverVelocityVerletPositionStep();
        inline void noseHooverVelocityVerletHalfVelocityStep(double& ke_half_sum, double& ke_sum);
        inline void noseHooverVelocityVerletFullVelocityStep();
        inline void setDpmPosition();
        inline void setDpmVelocity();
        inline void centerVelocities(double vx, double vy);
        inline void calcKineticEnergies();
        inline double getTemperature();
        inline double rescaleVelocitiesToTemp(double T_target);
        inline void setPos_i(const int i);
        inline void setPos_j(const int j);
        inline void setPos_k(const int k);
        inline void setProjLengthToPos_k(const double proj_length, const double bond_length);
        inline int getNextVertex(const int i);
        inline int getPrevVertex(const int i);
        inline std::vector<double> getVertexPos(const int i);

        // Member Functions
        void setAreaForceEnergy(const int i, const int k);
        void setBondBendStretchForcesEnergies(const int i, const int j, const int k);
        void innerDpmForceRoutine();
        void calcBondLengthsAnglesAreaPerimeter();
        void setParticleVertexForceEnergy(DPM2D& other_dpm, int other_vertex_i);
        void setParticleSegmentForceEnergy(DPM2D& other_dpm, int other_vertex_i);
        void setParticleVertexForceEnergyWCA(DPM2D& other_dpm, int other_vertex_i);
        void setParticleSegmentForceEnergyWCA(DPM2D& other_dpm, int other_vertex_i);
        void setInteractionForceEnergy(DPM2D& other_dpm);


        // Default Constructor
        DPM2D(double cx, double cy, double radius, double n_vertices, SimParams2D& simparams, double length_diam_ratio=2.0, double vx=0.0, double vy=0.0, double T_0=0.0, double seed=12345.0);
        
        // Default Destructor
        ~DPM2D();
};


// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- Inline Member Functions ------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

inline void DPM2D::verletPositionStep() {
    // update the positions
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->pos_vertex[this->simparams.N_dim * i + dim] += this->vel_vertex[this->simparams.N_dim * i + dim] * this->simparams.dt + 0.5 * this->force_vertex[this->simparams.N_dim * i + dim] * this->simparams.dt * this->simparams.dt / this->simparams.mass_vertex;
        }
    }
}

inline void DPM2D::verletVelocityStep() {
    // update the velocities and accelerations
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->vel_vertex[this->simparams.N_dim * i + dim] += 0.5 * this->simparams.dt * (this->force_vertex[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex + this->acc_vertex[this->simparams.N_dim * i + dim]);
            this->acc_vertex[this->simparams.N_dim * i + dim] = this->force_vertex[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex;
        }
    }
}

inline void DPM2D::noseHooverVelocityVerletPositionStep() {
    // update the positions
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->pos_vertex[this->simparams.N_dim * i + dim] += this->vel_vertex[this->simparams.N_dim * i + dim] * this->simparams.dt + 0.5 * (this->force_vertex[this->simparams.N_dim * i + dim] / simparams.mass_vertex - this->simparams.eta * this->vel_vertex[this->simparams.N_dim * i + dim]) * this->simparams.dt * this->simparams.dt;
        }
    }
}

inline void DPM2D::noseHooverVelocityVerletHalfVelocityStep(double& ke_half_sum, double& ke_sum) {
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            // store the half-step velocity in the acceleration for now
            this->acc_vertex[this->simparams.N_dim * i + dim] = this->vel_vertex[this->simparams.N_dim * i + dim] + 0.5 * (this->force_vertex[this->simparams.N_dim * i + dim] / simparams.mass_vertex - this->simparams.eta * this->vel_vertex[this->simparams.N_dim * i + dim]) * this->simparams.dt;
            // calculate the kinetic energy for the eta update
            ke_sum += this->vel_vertex[this->simparams.N_dim * i + dim] * this->vel_vertex[this->simparams.N_dim * i + dim] * this->simparams.mass_vertex / 2;
            // calculate the kinetic energy at the half-step for the upcoming eta update
            ke_half_sum += this->acc_vertex[this->simparams.N_dim * i + dim] * this->acc_vertex[this->simparams.N_dim * i + dim] * this->simparams.mass_vertex / 2;
        }
    }
}

inline void DPM2D::noseHooverVelocityVerletFullVelocityStep() {
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            // update the velocity using the half-step velocity (which is stored as the acceleration)
            this->vel_vertex[this->simparams.N_dim * i + dim] = (this->acc_vertex[this->simparams.N_dim * i + dim] + 0.5 * this->force_vertex[this->simparams.N_dim * i + dim] * this->simparams.dt / this->simparams.mass_vertex) / (1 + this->simparams.eta * this->simparams.dt / 2);
            // reset the acceleration using the full-step force (removes the half-step velocity from temporary storage here)
            this->acc_vertex[this->simparams.N_dim * i + dim] = this->force_vertex[this->simparams.N_dim * i + dim] / this->simparams.mass_vertex;
        }
    }
}

inline void DPM2D::setDpmPosition() {
    // sum over vertices and calculate the average position
    std::fill(this->pos_dpm.begin(), this->pos_dpm.end(), 0.0);
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->pos_dpm[dim] += this->pos_vertex[this->simparams.N_dim * i + dim];
        }
    }
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_dpm[dim] /= this->n_vertices;
    }
}

inline void DPM2D::setDpmVelocity() {
    // sum over vertices and calculate the average velocity
    std::fill(this->vel_dpm.begin(), this->vel_dpm.end(), 0.0);
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->vel_dpm[dim] += this->vel_vertex[this->simparams.N_dim * i + dim];
        }
    }
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->vel_dpm[dim] /= this->n_vertices;
    }
}

inline void DPM2D::centerVelocities(double vx, double vy) {
    this->setDpmVelocity();
    // subtract the average velocity from each velocity
    for (int i = 0; i < this->n_vertices; ++i) {
        this->vel_vertex[this->simparams.N_dim * i] -= this->vel_dpm[0] - vx;
        this->vel_vertex[this->simparams.N_dim * i + 1] -= this->vel_dpm[1] - vy;
    }
    // recalculate the dpm velocity
    this->setDpmVelocity();
}

inline void DPM2D::calcKineticEnergies() {
    // calculate the kinetic energy of the vertices
    this->kin_eng = 0.0;
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->kin_eng += this->simparams.mass_vertex * this->vel_vertex[this->simparams.N_dim * i + dim] * this->vel_vertex[this->simparams.N_dim * i + dim] / 2.0;
        }
    }
}

inline double DPM2D::getTemperature() {
    return 2.0 * this->kin_eng / (this->n_vertices * this->simparams.mass_vertex * this->simparams.kb);
}

inline double DPM2D::rescaleVelocitiesToTemp(double T_target) {
    this->calcKineticEnergies();
    double T_current = this->getTemperature();
    double scale_factor;
    if (T_current == 0.0) {
        scale_factor = 0.0;
    } else {
        scale_factor = std::sqrt(T_target / T_current);
    }
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
            this->vel_vertex[this->simparams.N_dim * i + dim] *= scale_factor;
        }
    }
    this->setDpmVelocity();
    this->calcKineticEnergies();
    return this->getTemperature();
}

inline void DPM2D::setPos_i(const int i) {
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_i[dim] = this->pos_vertex[this->simparams.N_dim * i + dim];
    }
}

inline void DPM2D::setPos_j(const int j) {
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_j[dim] = this->pos_vertex[this->simparams.N_dim * j + dim];
    }
}

inline void DPM2D::setPos_k(const int k) {
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_k[dim] = this->pos_vertex[this->simparams.N_dim * k + dim];
    }
}

inline void DPM2D::setProjLengthToPos_k(const double proj_length, const double bond_length) {  // requires that pos_j is set
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        this->pos_k[dim] = this->pos_j[dim] + proj_length * this->l_ij[dim] / bond_length;
    }
}

inline int DPM2D::getNextVertex(const int i) {
    return (i + 1) % this->n_vertices;
}

inline int DPM2D::getPrevVertex(const int i) {
    return (i - 1 + this->n_vertices) % this->n_vertices;
}

inline std::vector<double> DPM2D::getVertexPos(const int i) {
    std::vector<double> pos_i(this->simparams.N_dim);
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        pos_i[dim] = this->pos_vertex[this->simparams.N_dim * i + dim];
    }
    return pos_i;
}

#endif // DPM_HPP