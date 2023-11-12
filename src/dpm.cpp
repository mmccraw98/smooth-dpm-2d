#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <sstream>
#include <filesystem>

#include "dpm.hpp"
#include "misc.hpp"

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- DPM Class --------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //

// Assumes:
// All DPMS will have the same vertex diameter - sigma

DPM2D::DPM2D(double cx, double cy, double radius, double n_vertices, SimParams2D& simparams, double length_diam_ratio, double vx, double vy, double T_0, double seed) : simparams(simparams) {
    this->simparams = simparams;
    this->n_vertices = n_vertices;

    // fill vertex_is_active with true
    std::fill(this->vertex_is_active.begin(), this->vertex_is_active.end(), true);

    // calculate the dpm mass
    this->simparams.mass_dpm = this->simparams.mass_vertex * this->n_vertices;

    // initialize the vectors
    this->initializeVectors();

    // get the coordinates of the vertices
    this->pos_vertex = getCircleCoords(cx, cy, radius, n_vertices);
    this->setDpmPosition();

    // set the force interaction terms from the current configuration:
    this->calcBondLengthsAnglesAreaPerimeter();

    this->A_0 = this->area;
    this->l_0 = this->perimeter / this->n_vertices;
    this->theta_0 = 2 * M_PI / this->n_vertices;


    // set the sigma
    this->sigma = this->l_0 / length_diam_ratio;

    // initialize the velocities
    // ensures proper standard deviation
    std::vector<double> randomized_velocities = generateRandomNormalVector(this->n_vertices * this->simparams.N_dim, 0.0, 1.0, seed);
    this->centerVelocities(0.0, 0.0);
    double Temp = this->rescaleVelocitiesToTemp(T_0);

    // ensures proper average velocity
    this->centerVelocities(vx, vy);

    // set the forces and energies
    this->innerDpmForceRoutine();
}

DPM2D::~DPM2D() {
    // destructor    
}


void DPM2D::initializeVectors() {
    // resize all vectors to the correct size
    this->pos_vertex.resize(this->simparams.N_dim * this->n_vertices);
    this->vel_vertex.resize(this->simparams.N_dim * this->n_vertices);
    this->force_vertex.resize(this->simparams.N_dim * this->n_vertices);
    this->acc_vertex.resize(this->simparams.N_dim * this->n_vertices);
    this->bond_lengths.resize(this->n_vertices);
    this->bond_angles.resize(this->n_vertices);
    this->pos_dpm.resize(this->simparams.N_dim);
    this->vel_dpm.resize(this->simparams.N_dim);
    this->force_dpm.resize(this->simparams.N_dim);
    this->pos_i.resize(this->simparams.N_dim);
    this->pos_j.resize(this->simparams.N_dim);
    this->pos_k.resize(this->simparams.N_dim);
    this->l_ij.resize(this->simparams.N_dim);
    this->l_ki.resize(this->simparams.N_dim);
    this->vertex_is_active.resize(this->n_vertices);
}

// ------------------------------------------------------------------------------------------------------------------------------------------ //
// ---------------------------------------------------------- DPM Force Calculation --------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------------ //


void DPM2D::setAreaForceEnergy(const int i, const int k) {
    // calculate the magnitude of the interaction (multiply with constant for force / square and multiply for energy)
    double magnitude_area = (this->area - this->A_0);

    // distribute the forces
    // use the derivative of the trapezoid area formula with respect to pos_i to get force_i
    this->force_vertex[this->simparams.N_dim * i] -= this->simparams.k_a * magnitude_area * (this->pos_i[1] + this->pos_k[1]) / 2;
    this->force_vertex[this->simparams.N_dim * i + 1] -= this->simparams.k_a * magnitude_area * (this->pos_i[0] - this->pos_k[0]) / 2;
    // use the derivative of the trapezoid area formula with respect to pos_k to get force_k
    this->force_vertex[this->simparams.N_dim * k] += this->simparams.k_a * magnitude_area * (this->pos_i[1] + this->pos_k[1]) / 2;
    this->force_vertex[this->simparams.N_dim * k + 1] -= this->simparams.k_a * magnitude_area * (this->pos_i[0] - this->pos_k[0]) / 2;
    // NOTE - moved the potential energy update to the innerDpmForceRoutine() function
}

void DPM2D::setBondBendStretchForcesEnergies(const int i, const int j, const int k) {
    // calculate the bond length between i and j (update the bond_lengths vector at i)
    setDistVect(this->l_ij, this->pos_i, this->pos_j, this->simparams);
    this->bond_lengths[i] = getVectLength(this->l_ij, this->simparams);

    // calculate the bond length between i and k - this is doubly innefficient, but better than recalculating the distances for each force
    setDistVect(this->l_ki, this->pos_k, this->pos_i, this->simparams);
    this->bond_lengths[k] = getVectLength(this->l_ki, this->simparams);

    // calculate the dot products
    double C_ii = getDotProd(l_ij, l_ij, this->simparams);
    double C_kk = getDotProd(l_ki, l_ki, this->simparams);
    double C_ki = getDotProd(l_ki, l_ij, this->simparams);

    // calculate the angle between i, j, and k (update the bond_angles vector at i)
    this->bond_angles[i] = std::acos(C_ki / std::sqrt(C_ii * C_kk));

    // calculate the magnitude of each force (sqaure it to get energy)
    double magnitude_length = (this->bond_lengths[i] - this->l_0);
    double magnitude_angle = (this->bond_angles[i] - this->theta_0);
    double scaling_factor_angle = this->simparams.k_b * (this->bond_angles[i] - this->theta_0) / (std::sqrt(C_ii * C_kk) * std::sin(this->bond_angles[i]));

    // distribute the forces
    for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
        // length
        this->force_vertex[this->simparams.N_dim * i + dim] -= this->simparams.k_l * magnitude_length * l_ij[dim] / this->bond_lengths[i];
        this->force_vertex[this->simparams.N_dim * j + dim] += this->simparams.k_l * magnitude_length * l_ij[dim] / this->bond_lengths[i];

        // angle
        this->force_vertex[this->simparams.N_dim * i + dim] += scaling_factor_angle * (C_ki / C_kk * l_ki[dim] - C_ki / C_ii * l_ij[dim] + l_ki[dim] - l_ij[dim]);
        this->force_vertex[this->simparams.N_dim * j + dim] += scaling_factor_angle * (C_ki / C_ii * l_ij[dim] - l_ki[dim]);
        this->force_vertex[this->simparams.N_dim * k + dim] -= scaling_factor_angle * (C_ki / C_kk * l_ki[dim] - l_ij[dim]);
    }

    // increment the area and perimeter
    this->area += (this->pos_i[0] * this->pos_k[1] - this->pos_k[0] * this->pos_i[1]) / 2.0;
    this->perimeter += this->bond_lengths[i];

    // increment the energies
    this->pot_eng += magnitude_length * magnitude_length * this->simparams.k_l / 2.0 + magnitude_angle * magnitude_angle * this->simparams.k_b / 2.0;
}

void DPM2D::innerDpmForceRoutine() {
    // reset the forces and energies
    this->pot_eng = 0.0;
    this->kin_eng = 0.0;
    std::fill(this->force_dpm.begin(), this->force_dpm.end(), 0.0);
    std::fill(this->force_vertex.begin(), this->force_vertex.end(), 0.0);

    // reset the area and perimeter
    this->area = 0.0;
    this->perimeter = 0.0;

    // loop over all vertices to calculate the bending and stretching forces
    for (int i = 0; i < this->n_vertices; ++i) {
        // get the next and previous vertices
        int k = this->getNextVertex(i);
        int j = this->getPrevVertex(i);

        this->setPos_i(i);
        this->setPos_j(j);
        this->setPos_k(k);

        // set the forces and energies
        this->setBondBendStretchForcesEnergies(i, j, k);
        // NOTE cannot calculate the area here as it is being updated!
    }

    // loop over all vertices and calcualte the area force
    for (int i = 0; i < this->n_vertices; ++i) {
        // get the next and previous vertices
        int k = this->getNextVertex(i);

        this->setPos_i(i);
        this->setPos_k(k);

        // set the forces and energies
        this->setAreaForceEnergy(i, k);
    }
    double magnitude_area = (this->area - this->A_0);
    this->pot_eng += this->simparams.k_a * magnitude_area * magnitude_area / 2.0;  // NOTE - if this goes in the per-vertex force calculation, you need to divide by the number of vertices to avoid double counting
}

void DPM2D::calcBondLengthsAnglesAreaPerimeter() {
    // reset the area and perimeter
    this->area = 0.0;
    this->perimeter = 0.0;
    // loop over all vertices
    for (int i = 0; i < this->n_vertices; ++i) {
        // set the next and previous vertices
        int k = this->getNextVertex(i);
        int j = this->getPrevVertex(i);

        this->setPos_i(i);
        this->setPos_j(j);
        this->setPos_k(k);

        // calculate the bond length between i and j (update the bond_lengths vector at i)
        setDistVect(this->l_ij, this->pos_i, this->pos_j, this->simparams);
        this->bond_lengths[i] = getVectLength(this->l_ij, this->simparams);

        // calculate the bond length between i and k - this is doubly innefficient, but better than recalculating the distances for each force
        setDistVect(this->l_ki, this->pos_k, this->pos_i, this->simparams);
        this->bond_lengths[k] = getVectLength(this->l_ki, this->simparams);

        // calculate the dot products
        double C_ii = getDotProd(l_ij, l_ij, this->simparams);
        double C_kk = getDotProd(l_ki, l_ki, this->simparams);
        double C_ki = getDotProd(l_ki, l_ij, this->simparams);

        // calculate the angle between i, j, and k (update the bond_angles vector at i)
        this->bond_angles[i] = std::acos(C_ki / std::sqrt(C_ii * C_kk));

        // increment the area and perimeter
        this->area += (this->pos_i[0] * this->pos_k[1] - this->pos_k[0] * this->pos_i[1]) / 2.0;
        this->perimeter += this->bond_lengths[i];
    }
}

// TODO - unify the setParticleVertexForceEnergy and setParticleSegmentForceEnergy function variants to a single function which depends on the interaction type
// particle - particle force energy
void DPM2D::setParticleVertexForceEnergy(DPM2D& other_dpm, int other_vertex_i) {
    // set setPos_i(other_vertex_id) in the other dpm first!
    other_dpm.setPos_i(other_vertex_i);
    // loop over all vertices
    for (int i = 0; i < this->n_vertices; ++i) {
        if (this->vertex_is_active[i]) {
            // set this position
            this->setPos_i(i);
            // calculate the distance vector between the particle and the vertex
            setDistVect(this->l_ij, this->pos_i, other_dpm.pos_i, this->simparams);
            // calculate the distance between the particle and the vertex
            double dist = getVectLength(this->l_ij, this->simparams);
            // calculate the magnitude of the force
            double magnitude = (dist - this->sigma);
            if (dist < this->sigma) {  // if they touch, calculate the interaction
                // calculate the force
                for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
                    this->force_vertex[this->simparams.N_dim * i + dim] -= this->simparams.k * magnitude * this->l_ij[dim] / dist;
                    other_dpm.force_vertex[other_dpm.simparams.N_dim * other_vertex_i + dim] += this->simparams.k * magnitude * this->l_ij[dim] / dist;
                }
                // calculate the energy
                other_dpm.pot_eng += this->simparams.k * magnitude * magnitude / 2.0;
                this->pot_eng += this->simparams.k * magnitude * magnitude / 2.0;
            }
        }
        else {
            this->vertex_is_active[i] = true;
        }
    }
}

// particle - segment force energy
void DPM2D::setParticleSegmentForceEnergy(DPM2D& other_dpm, int other_vertex_i) {
    // set other positions
    other_dpm.setPos_i(other_vertex_i);

    // loop over all vertices
    for (int i = 0; i < this->n_vertices; ++i) {
        // get segment endpoint
        this->setPos_i(i);
        // get segment origin
        int j = this->getPrevVertex(i);
        this->setPos_j(j);
        // calculate the segment vector
        setDistVect(this->l_ij, this->pos_i, this->pos_j, this->simparams);
        // calculate the projection of the particle onto the segment
        double proj_length = calcProjection(other_dpm.pos_i, this->pos_j, this->l_ij, this->bond_lengths[i], this->simparams);
        // if the projection is within the segment, calculate the force
        if (proj_length > 0 and proj_length < this->bond_lengths[i]) {
            // calculate the position of the projection
            this->setProjLengthToPos_k(proj_length, bond_lengths[i]);
            // turn off the vertices
            this->vertex_is_active[i] = false;
            this->vertex_is_active[j] = false;

            // store the projection position in pos_k
            // calculate the distance vector between pos_k and other.pos_i
            setDistVect(this->l_ki, this->pos_k, other_dpm.pos_i, this->simparams);  // points to pos_k (projection on segment)

            // get the distance between the projection and the particle
            // calculate the distance between pos_k and other.pos_i
            double dist = getVectLength(this->l_ki, this->simparams);

            // calculate the overlap with the segment
            double magnitude = (dist - this->sigma);
            // TODO - this needs to be adjusted (and the vertex one too) so that it is positive when the particle is inside the DPM

            // if there is an overlap, calculate the force
            if (dist < this->sigma) {
                // the force is distributed between the two vertices using statics
                double ratio = proj_length / this->bond_lengths[i];
                for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
                    this->force_vertex[this->simparams.N_dim * i + dim] -= this->simparams.k * magnitude * ratio * this->l_ki[dim] / dist;
                    this->force_vertex[this->simparams.N_dim * j + dim] -= this->simparams.k * magnitude * (1.0 - ratio) * this->l_ki[dim] / dist;

                    other_dpm.force_vertex[other_dpm.simparams.N_dim * other_vertex_i + dim] += this->simparams.k * magnitude * this->l_ki[dim] / dist;
                }
                // calculate the energy
                other_dpm.pot_eng += this->simparams.k * magnitude * magnitude / 2.0;
                this->pot_eng += this->simparams.k * magnitude * magnitude / 2.0;
            }
        }
    }
}

// particle - particle force energy for WCA (repulsive lennard jones)
void DPM2D::setParticleVertexForceEnergyWCA(DPM2D& other_dpm, int other_vertex_i) {
    // set setPos_i(other_vertex_id) in the other dpm first!
    other_dpm.setPos_i(other_vertex_i);
    // loop over all vertices
    for (int i = 0; i < this->n_vertices; ++i) {
        if (this->vertex_is_active[i]) {
            // set this position
            this->setPos_i(i);
            // calculate the distance vector between the particle and the vertex
            setDistVect(this->l_ij, this->pos_i, other_dpm.pos_i, this->simparams);
            // calculate the distance between the particle and the vertex
            double dist = getVectLength(this->l_ij, this->simparams);

            if (dist < this->sigma * 1.122462048309373) {  // 2^(1/6) * sigma
                // calculate the force
                double sigma12 = std::pow(this->sigma, 12) / std::pow(dist, 12);
                double sigma6 = std::pow(this->sigma, 6) / std::pow(dist, 6);
                double force_magnitude = 24 * this->simparams.k * (2 * sigma12 - sigma6) / dist;
                double energy_magnitude = 4 * this->simparams.k * (sigma12 - sigma6) + this->simparams.k;

                for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
                    this->force_vertex[this->simparams.N_dim * i + dim] += force_magnitude * this->l_ij[dim] / dist;
                    other_dpm.force_vertex[other_dpm.simparams.N_dim * other_vertex_i + dim] -= force_magnitude * this->l_ij[dim] / dist;
                }
                // calculate the energy
                other_dpm.pot_eng += energy_magnitude / 2.0;
                this->pot_eng += energy_magnitude / 2.0;
            }
        }
        else {
            this->vertex_is_active[i] = true;
        }
    }
}

// particle - segment force energy for WCA (repulsive lennard jones)
void DPM2D::setParticleSegmentForceEnergyWCA(DPM2D& other_dpm, int other_vertex_i) {
    // set other positions
    other_dpm.setPos_i(other_vertex_i);

    // loop over all vertices
    for (int i = 0; i < this->n_vertices; ++i) {
        // get segment endpoint
        this->setPos_i(i);
        // get segment origin
        int j = this->getPrevVertex(i);
        this->setPos_j(j);
        // calculate the segment vector
        setDistVect(this->l_ij, this->pos_i, this->pos_j, this->simparams);
        // calculate the projection of the particle onto the segment
        double proj_length = calcProjection(other_dpm.pos_i, this->pos_j, this->l_ij, this->bond_lengths[i], this->simparams);
        // if the projection is within the segment, calculate the force
        if (proj_length > 0 and proj_length < this->bond_lengths[i]) {
            // calculate the position of the projection
            this->setProjLengthToPos_k(proj_length, bond_lengths[i]);
            // turn off the vertices
            this->vertex_is_active[i] = false;
            this->vertex_is_active[j] = false;

            // store the projection position in pos_k
            // calculate the distance vector between pos_k and other.pos_i
            setDistVect(this->l_ki, this->pos_k, other_dpm.pos_i, this->simparams);  // points to pos_k (projection on segment)

            // get the distance between the projection and the particle
            // calculate the distance between pos_k and other.pos_i
            double dist = getVectLength(this->l_ki, this->simparams);
            
            // if there is an overlap, calculate the force
            if (dist < this->sigma * 1.122462048309373) {  // 2^(1/6) * sigma

                double sigma12 = std::pow(this->sigma, 12) / std::pow(dist, 12);
                double sigma6 = std::pow(this->sigma, 6) / std::pow(dist, 6);
                double force_magnitude = 24 * this->simparams.k * (2 * sigma12 - sigma6) / dist;
                double energy_magnitude = 4 * this->simparams.k * (sigma12 - sigma6) + this->simparams.k;


                // the force is distributed between the two vertices using statics
                double ratio = proj_length / this->bond_lengths[i];
                for (int dim = 0; dim < this->simparams.N_dim; ++dim) {
                    this->force_vertex[this->simparams.N_dim * i + dim] += force_magnitude * ratio * this->l_ki[dim] / dist;  // MAY HAVE TO ADJUST SIGN
                    this->force_vertex[this->simparams.N_dim * j + dim] += force_magnitude * (1.0 - ratio) * this->l_ki[dim] / dist;  // MAY HAVE TO ADJUST SIGN

                    other_dpm.force_vertex[other_dpm.simparams.N_dim * other_vertex_i + dim] -= force_magnitude * this->l_ki[dim] / dist;  // MAY HAVE TO ADJUST SIGN
                }
                // calculate the energy  // DO THESE NEED TO BE DIVIDED BY 2?  yes
                other_dpm.pot_eng += energy_magnitude / 2.0;
                this->pot_eng += energy_magnitude / 2.0;
            }
        }
    }
}

// combine the two - NOTE - the inner DPM forces must have already been calculated
// first set the vertex_is_active all to true
// at the end, add the interaction potential to the total potential
// be careful to not double count the particle dpm interactions
void DPM2D::setInteractionForceEnergy(DPM2D& other_dpm) {
    std::fill(this->vertex_is_active.begin(), this->vertex_is_active.end(), true);  // this might be overly redundant but it fixes a bug in the tracer test
    // wherein if the tracer starts over a vertex, the forces arent calculated properly


    // loop over all the other dpm vertices
    for (int other_vertex_i = 0; other_vertex_i < other_dpm.n_vertices; ++other_vertex_i) {
        // calculate the particle - segment interaction
        this->setParticleSegmentForceEnergyWCA(other_dpm, other_vertex_i);
        // calculate the particle - vertex interaction
        this->setParticleVertexForceEnergyWCA(other_dpm, other_vertex_i);
    }
}
