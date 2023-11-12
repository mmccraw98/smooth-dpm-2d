#ifndef SIM_HPP
#define SIM_HPP

#include <cmath>
#include <vector>

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- Sim Params Struct ------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

struct SimParams2D {
    double box_size[2];  // Store x_length and y_length
    int N_dim;  // number of dimensions
    double dt;  // time step
    double kb;  // boltzmann constant
    double k;  // particle - particle interaction
    double k_l;  // bond stretching
    double k_b;  // bond bending
    double k_a;  // area constraint
    double mass_vertex;  // mass of the vertices
    double mass_dpm;  // mass of the dpm
    double sigma;  // particle - particle interaction distance (width of particles)
    double damping_coeff;  // damping coefficient
    double Q;  // for nose-hoover thermostat
    double eta; // for nose-hoover thermostat

    SimParams2D(double x_length, double y_length, double dt, double kb, double k, double k_l, double k_b, double k_a, double mass_vertex) {
        box_size[0] = x_length;
        box_size[1] = y_length;
        N_dim = 2;
        this->dt = dt;
        this->kb = kb;
        this->k = k;
        this->k_l = k_l;
        this->k_b = k_b;
        this->k_a = k_a;
        this->mass_vertex = mass_vertex;
        this->Q = 0.1;  // for nose-hoover thermostat
        this->eta = 1.0;  // for nose-hoover thermostat
        this->damping_coeff = 0.0;  // default to be 0!
        // these are default values and should be overwritten when the DPM2D is initiated
        double sigma = 1.0;  // particle - particle interaction distance (width of particles)
        double mass_dpm = 1.0;  // mass of the dpm
    }
};

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- Inline Distance Calculators --------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

inline double getDotProd(const std::vector<double>& vect1, const std::vector<double>& vect2, const SimParams2D& simparams) {
    double dot_prod = 0.0;
    for (int dim = 0; dim < simparams.N_dim; ++dim) {
        dot_prod += vect1[dim] * vect2[dim];
    };
    return dot_prod;
}

inline double getPbcDistPoint(const double& x1, const double& x2, const double& axis_length) {
    double dx = x1 - x2;
    dx -= axis_length * std::round(dx / axis_length);
    return dx;
}

inline double getDist(const double *this_pos, const double *other_pos, const SimParams2D& simparams) {
    double dist_sq = 0.0;
    double delta = 0.0;
    for (int dim = 0; dim < simparams.N_dim; ++dim) {
        delta = getPbcDistPoint(this_pos[dim], other_pos[dim], simparams.box_size[dim]);
        dist_sq += delta * delta;
    };
    return std::sqrt(dist_sq);
}

inline double getVectLength(const std::vector<double>& vect, const SimParams2D& simparams) {
    double vect_length_sq = 0.0;
    for (int dim = 0; dim < simparams.N_dim; ++dim) {
        vect_length_sq += vect[dim] * vect[dim];
    }
    return std::sqrt(vect_length_sq);
}

inline void setDistVect(std::vector<double>& dist_vect, const std::vector<double>& this_pos, const std::vector<double>& other_pos, const SimParams2D& simparams) {
    for (int dim = 0; dim < simparams.N_dim; ++dim) {
        dist_vect[dim] = getPbcDistPoint(this_pos[dim], other_pos[dim], static_cast<double>(simparams.box_size[dim]));
    }
}

inline std::vector<double> getDistVect(const std::vector<double>& this_pos, const std::vector<double>& other_pos, const SimParams2D& simparams) {
    std::vector<double> dist_vect(simparams.N_dim);
    for (int dim = 0; dim < simparams.N_dim; ++dim) {
        dist_vect[dim] = getPbcDistPoint(this_pos[dim], other_pos[dim], static_cast<double>(simparams.box_size[dim]));
    }
    return dist_vect;
}

inline double calcProjection(const std::vector<double> this_pos, const std::vector<double> segment_origin_pos, const std::vector<double> segment_vec, const double& segment_length, const SimParams2D& simparams) {
    return (getPbcDistPoint(this_pos[0], segment_origin_pos[0], simparams.box_size[0]) * segment_vec[0] + getPbcDistPoint(this_pos[1], segment_origin_pos[1], simparams.box_size[1]) * segment_vec[1]) / segment_length;
}

#endif // SIM_HPP