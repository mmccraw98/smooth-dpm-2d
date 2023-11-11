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

DPM2D::DPM2D(double cx, double cy, double radius, double n_vertices, GeomConfig2D& geomconfig, ForceParams forceparams, double length_diam_ratio, double vx, double vy, double T_0, double seed) : forceparams(forceparams), geomconfig(geomconfig) {
    this->forceparams = forceparams;
    this->geomconfig = geomconfig;
    this->n_vertices = n_vertices;

    // fill vertex_is_active with true
    std::fill(this->vertex_is_active.begin(), this->vertex_is_active.end(), true);

    // calculate the dpm mass
    this->forceparams.mass_dpm = this->forceparams.mass_vertex * this->n_vertices;

    // initialize the vectors
    this->initializeVectors();

    // get the coordinates of the vertices
    this->pos_vertex = getCircleCoords(cx, cy, radius, n_vertices);
    this->setDpmPosition();

    // set the force interaction terms from the current configuration:
    this->calcBondLengthsAnglesAreaPerimeter();

    this->forceparams.A_0 = this->area;
    this->forceparams.l_0 = this->perimeter / this->n_vertices;
    this->forceparams.theta_0 = 2 * M_PI / this->n_vertices;


    // set the sigma
    this->forceparams.sigma = this->forceparams.l_0 / length_diam_ratio;

    // initialize the velocities
    // TODO fix this so that we can assign a temperature (to the individual vertices) and a velocity to the dpm independently
    // ensures proper standard deviation
    this->assignRandomNormalVelocities(1.0, seed);
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
    this->pos_vertex.resize(this->geomconfig.N_dim * this->n_vertices);
    this->vel_vertex.resize(this->geomconfig.N_dim * this->n_vertices);
    this->force_vertex.resize(this->geomconfig.N_dim * this->n_vertices);
    this->acc_vertex.resize(this->geomconfig.N_dim * this->n_vertices);
    this->bond_lengths.resize(this->n_vertices);
    this->bond_angles.resize(this->n_vertices);
    this->pos_dpm.resize(this->geomconfig.N_dim);
    this->vel_dpm.resize(this->geomconfig.N_dim);
    this->force_dpm.resize(this->geomconfig.N_dim);
    this->pos_i.resize(this->geomconfig.N_dim);
    this->pos_j.resize(this->geomconfig.N_dim);
    this->pos_k.resize(this->geomconfig.N_dim);
    this->l_ij.resize(this->geomconfig.N_dim);
    this->l_ki.resize(this->geomconfig.N_dim);
    this->vertex_is_active.resize(this->n_vertices);
}

void DPM2D::setAreaForceEnergy(const int i, const int k) {
    // calculate the magnitude of the interaction (multiply with constant for force / square and multiply for energy)
    double magnitude_area = (this->area - this->forceparams.A_0);

    // distribute the forces
    // use the derivative of the trapezoid area formula with respect to pos_i to get force_i
    this->force_vertex[this->geomconfig.N_dim * i] -= this->forceparams.k_a * magnitude_area * (this->pos_i[1] + this->pos_k[1]) / 2;
    this->force_vertex[this->geomconfig.N_dim * i + 1] -= this->forceparams.k_a * magnitude_area * (this->pos_i[0] - this->pos_k[0]) / 2;
    // use the derivative of the trapezoid area formula with respect to pos_k to get force_k
    this->force_vertex[this->geomconfig.N_dim * k] += this->forceparams.k_a * magnitude_area * (this->pos_i[1] + this->pos_k[1]) / 2;
    this->force_vertex[this->geomconfig.N_dim * k + 1] -= this->forceparams.k_a * magnitude_area * (this->pos_i[0] - this->pos_k[0]) / 2;
    // NOTE - moved the potential energy update to the innerDpmForceRoutine() function
}

void DPM2D::setBondBendStretchForcesEnergies(const int i, const int j, const int k) {
    // calculate the bond length between i and j (update the bond_lengths vector at i)
    setDistVect(this->l_ij, this->pos_i, this->pos_j, this->geomconfig);
    this->bond_lengths[i] = getVectLength(this->l_ij, this->geomconfig);

    // calculate the bond length between i and k - this is doubly innefficient, but better than recalculating the distances for each force
    setDistVect(this->l_ki, this->pos_k, this->pos_i, this->geomconfig);
    this->bond_lengths[k] = getVectLength(this->l_ki, this->geomconfig);

    // calculate the dot products
    double C_ii = getDotProd(l_ij, l_ij, this->geomconfig);
    double C_kk = getDotProd(l_ki, l_ki, this->geomconfig);
    double C_ki = getDotProd(l_ki, l_ij, this->geomconfig);

    // calculate the angle between i, j, and k (update the bond_angles vector at i)
    this->bond_angles[i] = std::acos(C_ki / std::sqrt(C_ii * C_kk));

    // calculate the magnitude of each force (sqaure it to get energy)
    double magnitude_length = (this->bond_lengths[i] - this->forceparams.l_0);
    double magnitude_angle = (this->bond_angles[i] - this->forceparams.theta_0);
    double scaling_factor_angle = this->forceparams.k_b * (this->bond_angles[i] - this->forceparams.theta_0) / (std::sqrt(C_ii * C_kk) * std::sin(this->bond_angles[i]));

    // distribute the forces
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        // length
        this->force_vertex[this->geomconfig.N_dim * i + dim] -= this->forceparams.k_l * magnitude_length * l_ij[dim] / this->bond_lengths[i];
        this->force_vertex[this->geomconfig.N_dim * j + dim] += this->forceparams.k_l * magnitude_length * l_ij[dim] / this->bond_lengths[i];

        // angle
        this->force_vertex[this->geomconfig.N_dim * i + dim] += scaling_factor_angle * (C_ki / C_kk * l_ki[dim] - C_ki / C_ii * l_ij[dim] + l_ki[dim] - l_ij[dim]);
        this->force_vertex[this->geomconfig.N_dim * j + dim] += scaling_factor_angle * (C_ki / C_ii * l_ij[dim] - l_ki[dim]);
        this->force_vertex[this->geomconfig.N_dim * k + dim] -= scaling_factor_angle * (C_ki / C_kk * l_ki[dim] - l_ij[dim]);
    }

    // increment the area and perimeter
    this->area += (this->pos_i[0] * this->pos_k[1] - this->pos_k[0] * this->pos_i[1]) / 2.0;
    this->perimeter += this->bond_lengths[i];

    // increment the energies
    this->pot_eng += magnitude_length * magnitude_length * this->forceparams.k_l / 2.0 + magnitude_angle * magnitude_angle * this->forceparams.k_b / 2.0;
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
    double magnitude_area = (this->area - this->forceparams.A_0);
    this->pot_eng += this->forceparams.k_a * magnitude_area * magnitude_area / 2.0;  // NOTE - if this goes in the per-vertex force calculation, you need to divide by the number of vertices to avoid double counting
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
        setDistVect(this->l_ij, this->pos_i, this->pos_j, this->geomconfig);
        this->bond_lengths[i] = getVectLength(this->l_ij, this->geomconfig);

        // calculate the bond length between i and k - this is doubly innefficient, but better than recalculating the distances for each force
        setDistVect(this->l_ki, this->pos_k, this->pos_i, this->geomconfig);
        this->bond_lengths[k] = getVectLength(this->l_ki, this->geomconfig);

        // calculate the dot products
        double C_ii = getDotProd(l_ij, l_ij, this->geomconfig);
        double C_kk = getDotProd(l_ki, l_ki, this->geomconfig);
        double C_ki = getDotProd(l_ki, l_ij, this->geomconfig);

        // calculate the angle between i, j, and k (update the bond_angles vector at i)
        this->bond_angles[i] = std::acos(C_ki / std::sqrt(C_ii * C_kk));

        // increment the area and perimeter
        this->area += (this->pos_i[0] * this->pos_k[1] - this->pos_k[0] * this->pos_i[1]) / 2.0;
        this->perimeter += this->bond_lengths[i];
    }
}

void DPM2D::writeToLogFiles(std::ofstream& vertex_log, std::ofstream& dpm_log, int step, int dpm_id, int precision) {
    if (!vertex_log) {
        std::cout << "ERROR: could not write to vertex log file" << std::endl;
    }
    if (!dpm_log) {
        std::cout << "ERROR: could not write to dpm log file" << std::endl;
    }
    // fill force dpm, vel dpm, pos dpm with zeros
    std::fill(this->force_dpm.begin(), this->force_dpm.end(), 0.0);
    std::fill(this->vel_dpm.begin(), this->vel_dpm.end(), 0.0);
    std::fill(this->pos_dpm.begin(), this->pos_dpm.end(), 0.0);
    // set kinetic energies to zero
    this->kin_eng = 0.0;
    // loop over all vertices 
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            // add the vertex position to the dpm position
            this->pos_dpm[dim] += this->getVertexPos(i)[dim];
            // add the vertex velocity to the dpm velocity
            this->vel_dpm[dim] += this->vel_vertex[this->geomconfig.N_dim * i + dim];
            // add the vertex force to the dpm force
            this->force_dpm[dim] += this->force_vertex[this->geomconfig.N_dim * i + dim];
            // add the vertex kinetic energy to the dpm kinetic energy
            this->kin_eng += this->forceparams.mass_vertex * this->vel_vertex[this->geomconfig.N_dim * i + dim] * this->vel_vertex[this->geomconfig.N_dim * i + dim] / 2.0;
            // write the vertex log file
        }
        vertex_log << std::fixed << std::setprecision(precision) << step << "," << this->geomconfig.dt * step << "," << dpm_id << "," << i << "," << this->pos_vertex[this->geomconfig.N_dim * i] << "," << this->pos_vertex[this->geomconfig.N_dim * i + 1] << "," << this->vel_vertex[this->geomconfig.N_dim * i] << "," << this->vel_vertex[this->geomconfig.N_dim * i + 1] << "," << this->force_vertex[this->geomconfig.N_dim * i] << "," << this->force_vertex[this->geomconfig.N_dim * i + 1] << "\n";
    }

    // divide the dpm position by the number of vertices
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->pos_dpm[dim] /= this->n_vertices;
        this->vel_dpm[dim] /= this->n_vertices;
        this->force_dpm[dim] /= this->n_vertices;
    }
    // write the dpm log file
    dpm_log << std::fixed << std::setprecision(precision) << step << "," << this->geomconfig.dt * step << "," << dpm_id << "," << this->pos_dpm[0] << "," << this->pos_dpm[1] << "," << this->vel_dpm[0] << "," << this->vel_dpm[1] << "," << this->force_dpm[0] << "," << this->force_dpm[1] << "," << this->pot_eng << "," << this->kin_eng << "," << this->pot_eng + this->kin_eng << "," << this->area << "," << this->perimeter << "," << this->forceparams.mass_dpm << "," << this->forceparams.mass_vertex << "," << this->forceparams.sigma << "," << this->forceparams.l_0 << "," << this->forceparams.A_0 << "," << this->forceparams.theta_0 << "," << this->forceparams.k << "," << this->forceparams.k_l << "," << this->forceparams.k_b << "," << this->forceparams.k_a << "," << this->forceparams.damping_coeff << ","  << this->forceparams.Q << "," << this->forceparams.eta << "\n";
}

std::ofstream createVertexLogFile(const std::string& file_path) {
    std::ofstream file;
    file.open(file_path);
    if (!file.is_open()) {
        std::cout << "ERROR: could not open vertex log file " << file_path << std::endl;
    }
    file << "step,t,dpm_id,id,x,y,vx,vy,fx,fy\n";
    return file;
}

std::ofstream createDpmLogFile(const std::string& file_path) {
    std::ofstream file;
    file.open(file_path);
    if (!file.is_open()) {
        std::cout << "ERROR: could not open dpm log file " << file_path << std::endl;
    }
    file << "step,t,dpm_id,x,y,vx,vy,fx,fy,pe,ke,te,area,perim,mass,vertex_mass,vertex_sigma,length_0,area_0,angle_0,k,k_l,k_b,k_a,damping,Q,eta\n";
    return file;
}

void writeMacroConsoleHeader() {
    std::cout << std::string(12 * 8 + 6, '_') << std::endl;
    std::cout << std::setw(12) << "step" << " | " 
              << std::setw(12) << "time" << " | " 
              << std::setw(12) << "pe" << " | " 
              << std::setw(12) << "ke" << " | " 
              << std::setw(12) << "te" << " | "
              << std::setw(12) << "phi" << " | " 
              << std::setw(12) << "temp" << "\n";
    std::cout << std::string(12 * 8 + 6, '_') << std::endl;
}

std::ofstream createMacroLogFile(const std::string& file_path) {
    std::ofstream file;
    file.open(file_path);
    if (!file.is_open()) {
        std::cout << "ERROR: could not open macro log file " << file_path << std::endl;
    }
    file << "step,t,pe,ke,te,phi,temp,eta\n";

    return file;
}

std::ofstream createConfigLogFile(const std::string& file_path, GeomConfig2D& geomconfig) {
    // create the config file in the dir
    std::ofstream file;
    file.open(file_path);
    if (!file.is_open()) {
        std::cout << "ERROR: could not open macro log file " << file_path << std::endl;
    }
    file << "step,Lx,Ly,dt,kb" << std::endl;
    return file;
}

void writeToConfigLogFile(std::ofstream& config_log, GeomConfig2D& geomconfig, int step, int precision) {
    if (!config_log) {
        std::cout << "ERROR: could not write to config log file" << std::endl;
    }
    config_log << std::fixed << std::setprecision(precision) << step << "," << geomconfig.box_size[0] << "," << geomconfig.box_size[1] << "," << geomconfig.dt << "," << geomconfig.kb << "\n";
}

void writeToMacroLogFile(std::ofstream& macro_log, std::vector<DPM2D>& dpms, int step, GeomConfig2D& geomconfig, int precision, int console_log_freq) {
    if (!macro_log) {
        std::cout << "ERROR: could not write to macro log file" << std::endl;
    }
    // calculate the macro variables
    double pe = 0.0;
    double ke = 0.0;
    double phi = 0.0;
    int num_vertices = 0;
    for (int id = 0; id < dpms.size(); ++id) {
        pe += dpms[id].pot_eng;
        ke += dpms[id].kin_eng;
        phi += dpms[id].area;
        num_vertices += dpms[id].n_vertices;
    }
    // calculate the temperature
    double temp = 2 * ke / (geomconfig.N_dim * num_vertices * geomconfig.kb);
    // calculate the packing fraction
    phi /= (geomconfig.box_size[0] * geomconfig.box_size[1]);

    // write the macro log file
    macro_log << std::fixed << std::setprecision(precision) << step << "," << dpms[0].geomconfig.dt * step << "," << pe << "," << ke  << "," << pe + ke << "," << phi << "," << temp << "," << dpms[0].forceparams.eta << "\n";
    // log to console the same thing
    if (step % console_log_freq == 0) {
       std::cout << std::setw(12) << std::fixed << std::setprecision(3) << step << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << (dpms[0].geomconfig.dt * step) << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << pe << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << ke << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << ke + pe << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << phi << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << temp << "\n";
    }
}

std::vector<double> getCircleCoords(double cx, double cy, double radius, double n_vertices) {
    // reserve a vector of size 2*n_vertices and fill it with the x and y coordinates of the vertices
    std::vector<double> coords(2*n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        coords[2 * i] = cx + radius * cos(2 * M_PI * i / n_vertices);
        coords[2 * i + 1] = cy + radius * sin(2 * M_PI * i / n_vertices);
    }
    return coords;
}

void initDataFiles(std::string dir, GeomConfig2D geomconfig) {
    // make the dir
    std::string command = "mkdir -p " + dir;
    system(command.c_str());
}

void logDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step, int save_freq, int console_log_freq, std::ofstream& vertex_log, std::ofstream& dpm_log, std::ofstream& macro_log, std::ofstream& config_log, GeomConfig2D& geomconfig) {
    if (step % save_freq == 0) {
        for (int id = 0; id < num_dpms; ++id) {
            dpms[id].writeToLogFiles(vertex_log, dpm_log, step, id);
        }

        // write the macro variables to the log file
        writeToMacroLogFile(macro_log, dpms, step, geomconfig, 10, console_log_freq);
        writeToConfigLogFile(config_log, geomconfig, step, 10);
    }
}

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
            setDistVect(this->l_ij, this->pos_i, other_dpm.pos_i, this->geomconfig);
            // calculate the distance between the particle and the vertex
            double dist = getVectLength(this->l_ij, this->geomconfig);
            // calculate the magnitude of the force
            double magnitude = (dist - this->forceparams.sigma);
            if (dist < this->forceparams.sigma) {  // if they touch, calculate the interaction
                // calculate the force
                for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                    this->force_vertex[this->geomconfig.N_dim * i + dim] -= this->forceparams.k * magnitude * this->l_ij[dim] / dist;
                    other_dpm.force_vertex[other_dpm.geomconfig.N_dim * other_vertex_i + dim] += this->forceparams.k * magnitude * this->l_ij[dim] / dist;
                }
                // calculate the energy
                other_dpm.pot_eng += this->forceparams.k * magnitude * magnitude / 2.0;
                this->pot_eng += this->forceparams.k * magnitude * magnitude / 2.0;
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
        setDistVect(this->l_ij, this->pos_i, this->pos_j, this->geomconfig);
        // calculate the projection of the particle onto the segment
        double proj_length = calcProjection(other_dpm.pos_i, this->pos_j, this->l_ij, this->bond_lengths[i], this->geomconfig);
        // if the projection is within the segment, calculate the force
        if (proj_length > 0 and proj_length < this->bond_lengths[i]) {
            // calculate the position of the projection
            this->setProjLengthToPos_k(proj_length, bond_lengths[i]);
            // turn off the vertices
            this->vertex_is_active[i] = false;
            this->vertex_is_active[j] = false;

            // store the projection position in pos_k
            // calculate the distance vector between pos_k and other.pos_i
            setDistVect(this->l_ki, this->pos_k, other_dpm.pos_i, this->geomconfig);  // points to pos_k (projection on segment)

            // get the distance between the projection and the particle
            // calculate the distance between pos_k and other.pos_i
            double dist = getVectLength(this->l_ki, this->geomconfig);

            // calculate the overlap with the segment
            double magnitude = (dist - this->forceparams.sigma);
            // TODO - this needs to be adjusted (and the vertex one too) so that it is positive when the particle is inside the DPM

            // if there is an overlap, calculate the force
            if (dist < this->forceparams.sigma) {
                // the force is distributed between the two vertices using statics
                double ratio = proj_length / this->bond_lengths[i];
                for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                    this->force_vertex[this->geomconfig.N_dim * i + dim] -= this->forceparams.k * magnitude * ratio * this->l_ki[dim] / dist;
                    this->force_vertex[this->geomconfig.N_dim * j + dim] -= this->forceparams.k * magnitude * (1.0 - ratio) * this->l_ki[dim] / dist;

                    other_dpm.force_vertex[other_dpm.geomconfig.N_dim * other_vertex_i + dim] += this->forceparams.k * magnitude * this->l_ki[dim] / dist;
                }
                // calculate the energy
                other_dpm.pot_eng += this->forceparams.k * magnitude * magnitude / 2.0;
                this->pot_eng += this->forceparams.k * magnitude * magnitude / 2.0;
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
            setDistVect(this->l_ij, this->pos_i, other_dpm.pos_i, this->geomconfig);
            // calculate the distance between the particle and the vertex
            double dist = getVectLength(this->l_ij, this->geomconfig);

            if (dist < this->forceparams.sigma * 1.122462048309373) {  // 2^(1/6) * sigma
                // calculate the force
                double sigma12 = std::pow(this->forceparams.sigma, 12) / std::pow(dist, 12);
                double sigma6 = std::pow(this->forceparams.sigma, 6) / std::pow(dist, 6);
                double force_magnitude = 24 * this->forceparams.k * (2 * sigma12 - sigma6) / dist;
                double energy_magnitude = 4 * this->forceparams.k * (sigma12 - sigma6) + this->forceparams.k;

                for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                    this->force_vertex[this->geomconfig.N_dim * i + dim] += force_magnitude * this->l_ij[dim] / dist;
                    other_dpm.force_vertex[other_dpm.geomconfig.N_dim * other_vertex_i + dim] -= force_magnitude * this->l_ij[dim] / dist;
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
        setDistVect(this->l_ij, this->pos_i, this->pos_j, this->geomconfig);
        // calculate the projection of the particle onto the segment
        double proj_length = calcProjection(other_dpm.pos_i, this->pos_j, this->l_ij, this->bond_lengths[i], this->geomconfig);
        // if the projection is within the segment, calculate the force
        if (proj_length > 0 and proj_length < this->bond_lengths[i]) {
            // calculate the position of the projection
            this->setProjLengthToPos_k(proj_length, bond_lengths[i]);
            // turn off the vertices
            this->vertex_is_active[i] = false;
            this->vertex_is_active[j] = false;

            // store the projection position in pos_k
            // calculate the distance vector between pos_k and other.pos_i
            setDistVect(this->l_ki, this->pos_k, other_dpm.pos_i, this->geomconfig);  // points to pos_k (projection on segment)

            // get the distance between the projection and the particle
            // calculate the distance between pos_k and other.pos_i
            double dist = getVectLength(this->l_ki, this->geomconfig);
            
            // if there is an overlap, calculate the force
            if (dist < this->forceparams.sigma * 1.122462048309373) {  // 2^(1/6) * sigma

                double sigma12 = std::pow(this->forceparams.sigma, 12) / std::pow(dist, 12);
                double sigma6 = std::pow(this->forceparams.sigma, 6) / std::pow(dist, 6);
                double force_magnitude = 24 * this->forceparams.k * (2 * sigma12 - sigma6) / dist;
                double energy_magnitude = 4 * this->forceparams.k * (sigma12 - sigma6) + this->forceparams.k;


                // the force is distributed between the two vertices using statics
                double ratio = proj_length / this->bond_lengths[i];
                for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                    this->force_vertex[this->geomconfig.N_dim * i + dim] += force_magnitude * ratio * this->l_ki[dim] / dist;  // MAY HAVE TO ADJUST SIGN
                    this->force_vertex[this->geomconfig.N_dim * j + dim] += force_magnitude * (1.0 - ratio) * this->l_ki[dim] / dist;  // MAY HAVE TO ADJUST SIGN

                    other_dpm.force_vertex[other_dpm.geomconfig.N_dim * other_vertex_i + dim] -= force_magnitude * this->l_ki[dim] / dist;  // MAY HAVE TO ADJUST SIGN
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

void dampDpms(std::vector<DPM2D>& dpms, double damping) {
    for (int id = 0; id < dpms.size(); ++id) {
        for (int dim = 0; dim < dpms[id].geomconfig.N_dim; ++dim) {
            for (int i = 0; i < dpms[id].n_vertices; ++i) {
                dpms[id].force_vertex[dpms[id].geomconfig.N_dim * i + dim] -= damping * dpms[id].vel_vertex[dpms[id].geomconfig.N_dim * i + dim];
            }
        }
    }
}

void verletStepDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step, double damping) {
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].verletPositionStep();
    }

    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].innerDpmForceRoutine();
    }
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        for (int other_id = 0; other_id < num_dpms; ++other_id) {  // TODO - test if double counting
            if (id != other_id) {
                dpms[id].setInteractionForceEnergy(dpms[other_id]);
            }
        }
    }
    if (damping > 0) {
        dampDpms(dpms, damping);
    }

    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].verletVelocityStep();
    }
}

void noseHooverVelocityVerletStepDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step, double& eta, double T_target, double Q, double damping) {
    double eta_half = 0.0;
    double ke_half_sum = 0.0;
    double ke_sum = 0.0;
    
    // update the positions
    int num_particles = 0;
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].noseHooverVelocityVerletPositionStep(eta);
        num_particles += dpms[id].n_vertices;
    }
    double K = (dpms[0].geomconfig.N_dim * num_particles + 1) / 2 * dpms[0].geomconfig.kb * T_target;

    // update the half velocities and store them in acceleration temporarily (also updates the half-eta)
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].noseHooverVelocityVerletHalfVelocityStep(eta, Q, K, ke_half_sum, ke_sum);
    }

    // update the forces
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].innerDpmForceRoutine();
    }
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        for (int other_id = 0; other_id < num_dpms; ++other_id) {  // TODO - test if double counting
            if (id != other_id) {
                dpms[id].setInteractionForceEnergy(dpms[other_id]);
            }
        }
    }
    if (damping > 0) {
        dampDpms(dpms, damping);
    }
    
    // update the eta half-step using the kinetic energy sum
    eta_half = eta + dpms[0].geomconfig.dt / (2 * Q) * (ke_sum - K);

    // update the eta using the half-step kinetic energy sum
    eta = eta_half + dpms[0].geomconfig.dt / (2 * Q) * (ke_half_sum - K);

    // update the velocities and reset the acceleration using the force
    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].noseHooverVelocityVerletFullVelocityStep(eta);
    }
}


void tracerTest(std::string dir, int num_vertices) {
    GeomConfig2D geomconfig = GeomConfig2D(10.0, 10.0, 1e-3, 1.0);
    ForceParams forceparams = ForceParams(1.0, 1.0, 1.0, 1.0, 1.0);

    initDataFiles(dir, geomconfig);

    // make the logs
    std::ofstream macro_log = createMacroLogFile(dir + "macro_log.csv");
    std::ofstream vertex_log = createVertexLogFile(dir + "vertex_log.csv");
    std::ofstream dpm_log = createDpmLogFile(dir + "dpm_log.csv");

    writeMacroConsoleHeader();

    int num_dpms = 1;
    int save_freq = 100;
    int console_log_freq = 1000;
    double T_0 = 0.0;
    double vx_0 = 0.0;
    double vy_0 = 0.0;
    double cx = 5.0;
    double cy = 5.0;
    double R = 2.0;

    // tracer test

    int N_steps = 1000;

    // make the dpm
    DPM2D dpm = DPM2D(cx, cy, R, num_vertices, geomconfig, forceparams, 2.0, vx_0, vy_0, T_0, 12345.0);
    
    // make the tracer
    double trace_radius = R + dpm.forceparams.sigma / 2.0;
    DPM2D tracer = DPM2D(cx, cy, trace_radius, 1, geomconfig, forceparams, 2.0, vx_0, vy_0, T_0, 12345.0);
    tracer.forceparams.sigma = dpm.forceparams.sigma;

    // make a vector with two dpms
    std::vector<DPM2D> dpms;
    dpms.push_back(dpm);
    dpms.push_back(tracer);

    // loop
    for (int i = 0; i < N_steps; ++i) {
        // define dpms[1] (tracer) position as a function of theta
        double theta = 2.0 * M_PI * (i) / N_steps;
        dpms[1].pos_vertex[0] = cx + trace_radius * cos(theta);
        dpms[1].pos_vertex[1] = cy + trace_radius * sin(theta);

        // calculate the inner dpm forces
        dpms[0].innerDpmForceRoutine();
        dpms[1].force_vertex[0] = 0.0;
        dpms[1].force_vertex[1] = 0.0;

        // calculate the interaction forces between the dpm and the tracer
        dpms[0].setInteractionForceEnergy(dpms[1]);

        dpms[0].writeToLogFiles(vertex_log, dpm_log, i, 0);
        dpms[1].writeToLogFiles(vertex_log, dpm_log, i, 1);

        // std::cout << dpms[1].force_vertex[0] << " | " << dpms[1].force_vertex[1] << std::endl;
    }

    // close the logs
    vertex_log.close();
    dpm_log.close();
}


std::vector<double> generateLatticeCoordinates(int N, double lx, double ly) {
    std::vector<double> latticeCoordinates(2 * N, 0.0);

    int n_rows = static_cast<int>(std::sqrt(static_cast<double>(N)));
    int n_columns = n_rows; // Assuming a square grid for simplicity.

    double dx = lx / (n_columns + 1); // Gap between particles along x-axis.
    double dy = ly / (n_rows + 1);    // Gap between particles along y-axis.

    int count = 0;

    for (int i = 1; i <= n_rows; ++i) {
        for (int j = 1; j <= n_columns; ++j) {
            if (count >= N) {
                break; // Stop if we have reached the desired number of particles.
            }

            // Calculate lattice coordinates.
            double lattice_x = dx * j;
            double lattice_y = dy * i;

            // Store the coordinates into the vector.
            latticeCoordinates[2 * count] = lattice_x;
            latticeCoordinates[2 * count + 1] = lattice_y;

            count++;
        }
    }

    return latticeCoordinates;
}


double setLinearHarmonicForces(std::vector<double>& force, std::vector<double>& pos, std::vector<double>& distance_vector, int num_disks, GeomConfig2D& geomconfig, ForceParams& forceparams, std::vector<double> radii_list) {
    // calculate the area of all the particles
    double area = 0.0;
    
    // define a force routine
    for (int i = 0; i < num_disks; ++i) {
        area += M_PI * radii_list[i] * radii_list[i];
        for (int j = i + 1; j < num_disks; ++j) {
            // calculate the distance between the two disks
            for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                distance_vector[dim] = getPbcDistPoint(pos[i * geomconfig.N_dim + dim], pos[j * geomconfig.N_dim + dim], geomconfig.box_size[dim]);
            }
            double length = getVectLength(distance_vector, geomconfig);
            double overlap = (radii_list[i] + radii_list[j]) - length;  // NOTE - in bidisperse case, this will be (sigma_i + sigma_j) / 2 - length
            if (overlap > 0) {
                // calculate the force
                for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                    double f = forceparams.k * overlap * distance_vector[dim] / length;
                    force[i * geomconfig.N_dim + dim] += f;
                    force[j * geomconfig.N_dim + dim] -= f;
                }

                // subtract the overlap from the area
                // NOTE will need to change this for bidisperse case
                // https://mathworld.wolfram.com/Circle-CircleIntersection.html
                // double R = forceparams.sigma / 2.0;
                // double overlap_area = 2 * (R * R * std::acos(overlap / (2 * R))) - (overlap / 2) * std::sqrt(4 * R * R - overlap * overlap);
                // area -= overlap_area;
            }
        }
    }
    return area;
}


double adamMinimizeDpmForce(std::vector<DPM2D> dpms, GeomConfig2D& geomconfig, int max_steps, double alpha, double beta1, double beta2, double epsilon) {
    // will store first and second moments for each dpm in their acceleration and velocity vectors
    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[i].n_vertices; ++i) {
            for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                dpms[id].acc_vertex[geomconfig.N_dim * i + dim] = 0.0;
                dpms[id].vel_vertex[geomconfig.N_dim * i + dim] = 0.0;
                dpms[id].force_vertex[geomconfig.N_dim * i + dim] = 0.0;
            }
        }
    }

    double area = 0.0;
    double force_norm = 0.0;

    for (int t = 0; t < max_steps; ++t) {
        // calculate the forces
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

        // update the moments
        for (int id = 0; id < dpms.size(); ++id) {
            for (int i = 0; i < dpms[i].n_vertices; ++i) {
                for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                    // update the first moment
                    // NOTE the negative sign on the gradient is because we get the gradient of the potential from the force (which is the negative gradient of the potential)
                    dpms[id].acc_vertex[geomconfig.N_dim * i + dim] = beta1 * dpms[id].acc_vertex[geomconfig.N_dim * i + dim] - (1 - beta1) * dpms[id].force_vertex[geomconfig.N_dim * i + dim];
                    
                    // update the second moment
                    // NOTE - no need for the negative sign here because we are squaring the gradient
                    dpms[id].vel_vertex[geomconfig.N_dim * i + dim] = beta2 * dpms[id].vel_vertex[geomconfig.N_dim * i + dim] + (1 - beta2) * dpms[id].force_vertex[geomconfig.N_dim * i + dim] * dpms[id].force_vertex[geomconfig.N_dim * i + dim];

                    // calculate the bias corrected first moment
                    double first_moment_bias_corrected = dpms[id].acc_vertex[geomconfig.N_dim * i + dim] / (1 - std::pow(beta1, t + 1));

                    // calculate the bias corrected second moment
                    double second_moment_bias_corrected = dpms[id].vel_vertex[geomconfig.N_dim * i + dim] / (1 - std::pow(beta2, t + 1));

                    // update the position
                    dpms[id].pos_vertex[geomconfig.N_dim * i + dim] -= alpha * first_moment_bias_corrected / (std::sqrt(second_moment_bias_corrected) + epsilon);

                    // update the force norm
                    force_norm += std::pow(dpms[id].force_vertex[geomconfig.N_dim * i + dim], 2.0);
                }
            }
        }

        // calculate the force norm
        force_norm = std::sqrt(force_norm);
        std::cout << t << " " << force_norm << std::endl;

        // check if the force norm is below the threshold
        if (force_norm < epsilon) {
            std::cout << force_norm << " success" << std::endl;
            break;
        }

        // reset all values to 0!
        for (int id = 0; id < dpms.size(); ++id) {
            for (int i = 0; i < dpms[i].n_vertices; ++i) {
                for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                    dpms[id].acc_vertex[geomconfig.N_dim * i + dim] = 0.0;
                    dpms[id].vel_vertex[geomconfig.N_dim * i + dim] = 0.0;
                    dpms[id].force_vertex[geomconfig.N_dim * i + dim] = 0.0;
                }
            }
        }
    }

    // calculate the area
    for (int id = 0; id < dpms.size(); ++id) {
        area += dpms[id].area;
    }
    return area / (geomconfig.box_size[0] * geomconfig.box_size[1]);  // TODO - generalize this
}

double adamMinimizeDiskForces(std::vector<double>& pos, GeomConfig2D& geomconfig, ForceParams& forceparams, std::vector<double>& radii_list, int num_disks, int max_steps, double alpha, double beta1, double beta2, double epsilon) {
    std::vector<double> gradient, first_moment, second_moment;  // these are to store the adam optimizer values
    std::vector<double> distance_vector(geomconfig.N_dim, 0.0);

    for (int i = 0; i < num_disks; ++i) {
        for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
            gradient.push_back(0.0);
            first_moment.push_back(0.0);
            second_moment.push_back(0.0);
        }
    }

    double force_norm = 0.0;
    double area = 0.0;

    // run the loop
    for (int t = 0; t < max_steps; ++t) {
        // calculate the forces
        area = setLinearHarmonicForces(gradient, pos, distance_vector, num_disks, geomconfig, forceparams, radii_list);
        // set the force norm back to 0
        force_norm = 0.0;

        // update the moments
        for (int i = 0; i < num_disks; ++i) {
            for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                // update the first moment
                // NOTE the negative sign on the gradient is because we get the gradient of the potential from the force (which is the negative gradient of the potential)
                first_moment[i * geomconfig.N_dim + dim] = beta1 * first_moment[i * geomconfig.N_dim + dim] - (1 - beta1) * gradient[i * geomconfig.N_dim + dim];
                
                // update the second moment
                // NOTE - no need for the negative sign here because we are squaring the gradient
                second_moment[i * geomconfig.N_dim + dim] = beta2 * second_moment[i * geomconfig.N_dim + dim] + (1 - beta2) * gradient[i * geomconfig.N_dim + dim] * gradient[i * geomconfig.N_dim + dim];

                // calculate the bias corrected first moment
                double first_moment_bias_corrected = first_moment[i * geomconfig.N_dim + dim] / (1 - std::pow(beta1, t + 1));

                // calculate the bias corrected second moment
                double second_moment_bias_corrected = second_moment[i * geomconfig.N_dim + dim] / (1 - std::pow(beta2, t + 1));

                // update the position
                pos[i * geomconfig.N_dim + dim] -= alpha * first_moment_bias_corrected / (std::sqrt(second_moment_bias_corrected) + epsilon);

                // enforce periodic boundary conditions NOTE might not be necessary
                if (pos[i * geomconfig.N_dim + dim] > geomconfig.box_size[dim]) {
                    pos[i * geomconfig.N_dim + dim] -= geomconfig.box_size[dim];
                } else if (pos[i * geomconfig.N_dim + dim] < 0.0) {
                    pos[i * geomconfig.N_dim + dim] += geomconfig.box_size[dim];
                }

                // update the force norm
                force_norm += std::pow(gradient[i * geomconfig.N_dim + dim], 2.0);
            }
        }

        // calculate the force norm
        force_norm = std::sqrt(force_norm);

        // check if the force norm is below the threshold
        if (force_norm < epsilon) {
            // std::cout << force_norm << " success" << std::endl;
            break;
        }

        // reset the gradient and moments
        std::fill_n(gradient.begin(), gradient.size(), 0.0);
        std::fill_n(first_moment.begin(), first_moment.size(), 0.0);
        std::fill_n(second_moment.begin(), second_moment.size(), 0.0);
    }
    return area / (geomconfig.box_size[0] * geomconfig.box_size[1]);  // TODO - generalize this
}




void DPM2D::innerDpmForceRoutinePlateCompression(std::vector<double> wall_bounds, double wall_strength, std::vector<double>& force_area, std::vector<double>& boundary_pos) {
    for (int i = 0; i < this->n_vertices; ++i) {
        if (this->pos_vertex[geomconfig.N_dim * i] < wall_bounds[0]) {
            this->force_vertex[geomconfig.N_dim * i] += wall_strength * (wall_bounds[0] - this->pos_vertex[geomconfig.N_dim * i]);
            force_area[0] += std::abs(wall_strength * (wall_bounds[0] - this->pos_vertex[geomconfig.N_dim * i]));
            boundary_pos[i] = this->pos_vertex[geomconfig.N_dim * i + 1];  // TODO make this more general - corners will ruin it
        }
        if (this->pos_vertex[geomconfig.N_dim * i] > wall_bounds[1]) {
            this->force_vertex[geomconfig.N_dim * i] += wall_strength * (wall_bounds[1] - this->pos_vertex[geomconfig.N_dim * i]);
            force_area[0] += std::abs(wall_strength * (wall_bounds[1] - this->pos_vertex[geomconfig.N_dim * i]));
            boundary_pos[i] = this->pos_vertex[geomconfig.N_dim * i + 1];  // TODO make this more general - corners will ruin it
        }
        if (this->pos_vertex[geomconfig.N_dim * i + 1] < wall_bounds[2]) {
            this->force_vertex[geomconfig.N_dim * i + 1] += wall_strength * (wall_bounds[2] - this->pos_vertex[geomconfig.N_dim * i + 1]);
            force_area[0] += std::abs(wall_strength * (wall_bounds[2] - this->pos_vertex[geomconfig.N_dim * i + 1]));
            boundary_pos[i] = this->pos_vertex[geomconfig.N_dim * i + 1];  // TODO make this more general - corners will ruin it
        }
        if (this->pos_vertex[geomconfig.N_dim * i + 1] > wall_bounds[3]) {
            this->force_vertex[geomconfig.N_dim * i + 1] += wall_strength * (wall_bounds[3] - this->pos_vertex[geomconfig.N_dim * i + 1]);
            force_area[0] += std::abs(wall_strength * (wall_bounds[3] - this->pos_vertex[geomconfig.N_dim * i + 1]));
            boundary_pos[i] = this->pos_vertex[geomconfig.N_dim * i + 1];  // TODO make this more general - corners will ruin it
        }
    }
}

void DPM2D::gradDescMinDpm(int max_steps, double eta, double tol, std::vector<double> wall_pos) {
    std::fill(this->force_vertex.begin(), this->force_vertex.end(), 0.0);
    std::vector<bool> vertex_is_fixed(this->n_vertices, false);
    double force_norm_prev = 0.0;
    double force_norm_curr = 0.0;

    for (int t = 0; t < max_steps; ++t) {
        // calculate the forces
        this->innerDpmForceRoutine();

        // apply fixes to the vertices
        for (int i = 0; i < this->n_vertices; ++i) {
            if (this->pos_vertex[geomconfig.N_dim * i] <= wall_pos[0]) {
                this->force_vertex[geomconfig.N_dim * i] += 1000.0 * std::pow((wall_pos[0] - this->pos_vertex[geomconfig.N_dim * i]), 2);
                // this->pos_vertex[geomconfig.N_dim * i] = wall_pos[0];
                // vertex_is_fixed[i] = true;
            }
            if (this->pos_vertex[geomconfig.N_dim * i] >= wall_pos[1]) {
                this->force_vertex[geomconfig.N_dim * i] -= 1000.0 * std::pow((wall_pos[1] - this->pos_vertex[geomconfig.N_dim * i]), 2);
                // this->pos_vertex[geomconfig.N_dim * i] = wall_pos[1];
                // vertex_is_fixed[i] = true;
            }
        }

        // apply the position updates to the non-fixed vertices
        force_norm_curr = 0.0;
        for (int i =0; i < this->n_vertices; ++i) {
            for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                // if (not vertex_is_fixed[i] or dim > 0) {
                double pos_update = eta * this->force_vertex[i * this->geomconfig.N_dim + dim];
                this->pos_vertex[i * this->geomconfig.N_dim + dim] += pos_update;
                // }
                force_norm_curr += this->force_vertex[i * this->geomconfig.N_dim + dim] * this->force_vertex[i * this->geomconfig.N_dim + dim];
            }
        }

        // if the force norm is nan break
        if (std::isnan(force_norm_curr)) {
            std::cout << "nan" << std::endl;
            return;
        }

        // check if the difference between the force norms is below the tolerance
        if (std::abs(force_norm_curr - force_norm_prev) / force_norm_curr < tol) {
            // std::cout << force_norm_curr - force_norm_prev << " success" << std::endl;
            return;
        }

        force_norm_prev = force_norm_curr;
    }
    std::cout << " failed to converge" << std::endl;
}

void DPM2D::gradDescMinDpmPlateForces(int max_steps, double eta, double tol, std::vector<double> wall_bounds, double wall_strength, std::vector<double>& force_area, std::vector<double>& boundary_pos) {
    std::fill(this->force_vertex.begin(), this->force_vertex.end(), 0.0);
    std::fill(force_area.begin(), force_area.end(), 0.0);
    std::fill(boundary_pos.begin(), boundary_pos.end(), 0.0);
    
    double disp_norm = 0.0;

    for (int t = 0; t < max_steps; ++t) {
        this->innerDpmForceRoutine();
        this->innerDpmForceRoutinePlateCompression(wall_bounds, wall_strength, force_area, boundary_pos);

        disp_norm = 0.0;
        for (int i = 0; i < this->n_vertices; ++i) {
            for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                double pos_update = eta * this->force_vertex[i * this->geomconfig.N_dim + dim];
                this->pos_vertex[i * this->geomconfig.N_dim + dim] += pos_update;
                disp_norm += pos_update * pos_update;
            }
        }

        if (disp_norm < tol) {
            // std::cout << disp_norm << " success" << std::endl;
            std::fill(force_area.begin(), force_area.end(), 0.0);
            std::fill(boundary_pos.begin(), boundary_pos.end(), 0.0);
            this->innerDpmForceRoutinePlateCompression(wall_bounds, wall_strength, force_area, boundary_pos);
            return;
        }
    }
    std::fill(force_area.begin(), force_area.end(), 0.0);
    std::fill(boundary_pos.begin(), boundary_pos.end(), 0.0);
    this->innerDpmForceRoutinePlateCompression(wall_bounds, wall_strength, force_area, boundary_pos);
            
    std::cout << disp_norm << " failed to converge" << std::endl;
}


void DPM2D::adamMinimizeDpmPlateForces(int max_steps, double alpha, double beta1, double beta2, double epsilon, std::vector<double> wall_bounds, double wall_strength) {
    // storing the gradient in force
    std::fill(this->force_vertex.begin(), this->force_vertex.end(), 0.0);
    // storing the first moment in acceleration
    std::fill(this->acc_vertex.begin(), this->acc_vertex.end(), 0.0);
    // storing the second moment in velocity
    std::fill(this->vel_vertex.begin(), this->vel_vertex.end(), 0.0);

    double force_norm = 0.0;

    std::vector<double> force_area(1, 0.0);
    std::vector<double> boundary_pos(this->n_vertices, 0.0);

    // run the loop
    for (int t = 0; t < max_steps; ++t) {
        // calculate the forces
        this->innerDpmForceRoutine();
        this->innerDpmForceRoutinePlateCompression(wall_bounds, wall_strength, force_area, boundary_pos);

        // set the force norm back to 0
        force_norm = 0.0;

        // update the moments
        for (int i = 0; i < this->n_vertices; ++i) {
            for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
                // update the first moment
                // NOTE the negative sign on the gradient is because we get the gradient of the potential from the force (which is the negative gradient of the potential)
                this->acc_vertex[i * this->geomconfig.N_dim + dim] = beta1 * this->acc_vertex[i * this->geomconfig.N_dim + dim] - (1 - beta1) * this->force_vertex[i * this->geomconfig.N_dim + dim];
                
                // update the second moment
                // NOTE - no need for the negative sign here because we are squaring the gradient
                this->vel_vertex[i * this->geomconfig.N_dim + dim] = beta2 * this->vel_vertex[i * this->geomconfig.N_dim + dim] + (1 - beta2) * this->force_vertex[i * this->geomconfig.N_dim + dim] * this->force_vertex[i * this->geomconfig.N_dim + dim];

                // calculate the bias corrected first moment
                double first_moment_bias_corrected = this->acc_vertex[i * this->geomconfig.N_dim + dim] / (1 - std::pow(beta1, t + 1));

                // calculate the bias corrected second moment
                double second_moment_bias_corrected = this->vel_vertex[i * this->geomconfig.N_dim + dim] / (1 - std::pow(beta2, t + 1));

                // update the position
                this->pos_vertex[i * this->geomconfig.N_dim + dim] -= alpha * first_moment_bias_corrected / (std::sqrt(second_moment_bias_corrected) + epsilon);

                // enforce periodic boundary conditions NOTE might not be necessary
                // if (this->pos_vertex[i * this->geomconfig.N_dim + dim] > this->geomconfig.box_size[dim]) {
                //     this->pos_vertex[i * this->geomconfig.N_dim + dim] -= this->geomconfig.box_size[dim];
                // } else if (this->pos_vertex[i * this->geomconfig.N_dim + dim] < 0.0) {
                //     this->pos_vertex[i * this->geomconfig.N_dim + dim] += this->geomconfig.box_size[dim];
                // }

                // update the force norm
                force_norm += std::pow(this->force_vertex[i * this->geomconfig.N_dim + dim], 2.0);
            }
        }

        // calculate the force norm
        force_norm = std::sqrt(force_norm);


        // check if the force norm is below the threshold
        if (force_norm < epsilon) {
            this->kin_eng = force_norm;
            std::cout << t << std::endl;
            return;
        }

        // reset the gradient and moments
        std::fill_n(this->force_vertex.begin(), this->force_vertex.size(), 0.0);
        std::fill_n(this->acc_vertex.begin(), this->acc_vertex.size(), 0.0);
        std::fill_n(this->vel_vertex.begin(), this->vel_vertex.size(), 0.0);
    }
    std::cout << "failed to converge " << force_norm << " " << epsilon << std::endl;
}

void plateCompressionSweep(std::string dir_base, int num_vertices, int N_points) {
    // define the dir where everything will be saved
    int iter = 0;

    double lower_bound = std::log(1);
    double upper_bound = std::log(500);
    double step_size = (upper_bound - lower_bound) / (N_points - 1);

    for (int i = 0; i < N_points; ++i) {
        for (int j = 0; j < N_points; ++j) {
            for (int k = 0; k < N_points; ++k) {

                // single dpm parallel plate compression test
                iter += 1;

                double k_a = std::floor(std::exp(lower_bound + i * step_size));
                double k_l = std::floor(std::exp(lower_bound + j * step_size));
                double k_b = std::floor(std::exp(lower_bound + k * step_size));

                GeomConfig2D geomconfig = GeomConfig2D(10.0, 10.0, 1e-3, 1.0);
                ForceParams forceparams = ForceParams(200.0, k_l, k_b, k_a, 1.0);

                std::string dir = dir_base + std::to_string(iter) + "/";
                std::cout << dir << " of " << std::pow(N_points, 3.0) << std::endl;

                initDataFiles(dir, geomconfig);

                // make the logs
                std::ofstream macro_log = createMacroLogFile(dir + "macro_log.csv");
                std::ofstream vertex_log = createVertexLogFile(dir + "vertex_log.csv");
                std::ofstream dpm_log = createDpmLogFile(dir + "dpm_log.csv");
                std::ofstream config_log = createConfigLogFile(dir + "config_log.csv", geomconfig);

                writeMacroConsoleHeader();

                int num_steps = 20000;
                int save_freq = 10;
                int console_log_freq = 500;
                double R = 2.0;
                double dr = 0.00001;

                DPM2D dpm = DPM2D(5.0, 5.0, R, num_vertices, geomconfig, forceparams, 2.0, 0.0, 0.0, 0.0, 0.0);
                std::vector<DPM2D> dpms;
                dpms.push_back(dpm);

                std::vector<double> wall_pos = {2.99, 7.01};
                // relax the dpm within the wall boundaries

                for (int step = 0; step < num_steps; ++step) {
                    dpms[0].gradDescMinDpm(200000, 0.0000001, 1e-5, wall_pos);
                    if (step % save_freq == 0) {
                        logDpmList(dpms.size(), dpms, step, save_freq, console_log_freq, vertex_log, dpm_log, macro_log, config_log, geomconfig);
                    }
                    wall_pos[0] += dr;
                    wall_pos[1] -= dr;
                }

                vertex_log.close();
                dpm_log.close();
                config_log.close();
                macro_log.close();
            }
        }
    }
}

PosRadius generateDiskPackCoords(int num_disks, std::vector<double> radii, std::vector<double> fraction, GeomConfig2D& geomconfig, ForceParams& forceparams, double seed, double dr, double tol, double phi_target, double num_steps) {
    double sum = 0.0;
    int total = 0;
    for (int i = 0; i < fraction.size(); i++) {
        sum += fraction[i];
        total += 1;
    }
    if (sum != 1.0) {
        std::cout << "ERROR: fraction does not add to 1.0" << std::endl;
    }

    // make the radius list
    std::vector<double> radii_list;
    for (int i = 0; i < fraction.size(); i++) {
        int num = std::round(fraction[i] * num_disks);
        for (int j = 0; j < num; j++) {
            radii_list.push_back(radii[i]);
        }
    }
    num_disks = radii_list.size();
    std::cout << "Made " << num_disks << " disks" << std::endl;

    // initialize the positions of the dpms using disks to prevent overlaps
    std::vector<double> pos;
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < radii_list.size(); ++i) {
        for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
            pos.push_back(dis(gen) * geomconfig.box_size[dim]);
        }
    }
    for (int step = 0; step < num_steps; ++step) {
        double phi = adamMinimizeDiskForces(pos, geomconfig, forceparams, radii_list, num_disks, 100000, 0.01, 0.9, 0.999, tol);
        geomconfig.box_size[0] -= dr;
        geomconfig.box_size[1] -= dr;
        if (step % 100 == 0) {
            std::cout << geomconfig.box_size[0] << " " << phi << std::endl;
        }
        if (phi > phi_target) {
            std::cout << geomconfig.box_size[0] << " " << phi << std::endl;
            break;
        }
    }
    // make the pos_radius struct
    PosRadius pos_rad;
    pos_rad.pos = pos;
    pos_rad.radii = radii_list;
    return pos_rad;
}

std::vector<DPM2D> generateDpmsFromDiskPack(PosRadius& pos_rad, GeomConfig2D& geomconfig, ForceParams& forceparams, double vertex_circumferencial_density, double radius_shrink_factor) {
    // assign the coordinates to the dpms
    std::vector<DPM2D> dpms;
    for (int i = 0; i < pos_rad.radii.size(); ++i) {
        double radius = pos_rad.radii[i] * radius_shrink_factor;  // to be really sure there are no overlaps
        int num_vertices = std::round(vertex_circumferencial_density * 2 * M_PI * radius);
        dpms.push_back(DPM2D(pos_rad.pos[2*i], pos_rad.pos[2*i+1], radius, num_vertices, geomconfig, forceparams, 2.0, 0.0, 0.0, 0.0, 0.0));
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

    // for (int id = 0; id < dpms.size(); ++id) {
    //     double vx_mean = 0.0;
    //     double vy_mean = 0.0;
        
    //     for (int i = 0; i < dpms[id].n_vertices; ++i) {
    //         vx_mean += dpms[id].vel_vertex[2 * i];
    //         vy_mean += dpms[id].vel_vertex[2 * i + 1];
    //     }

    //     vx_mean /= dpms[id].n_vertices;
    //     vy_mean /= dpms[id].n_vertices;

    //     for (int i = 0; i < dpms[id].n_vertices; ++i) {
    //         dpms[id].vel_vertex[2 * i] += vx - vx_mean;
    //         dpms[id].vel_vertex[2 * i + 1] += vy - vy_mean;
    //     }
    // }
}

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

void scaleDpmsToTemp(std::vector<DPM2D>& dpms, GeomConfig2D& geomconfig, ForceParams& forceparams, double temp_target, double seed) {
    // randomly assign velocities
    srand48(seed);
    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            dpms[id].vel_vertex[2 * i] = (drand48() - 0.5);
            dpms[id].vel_vertex[2 * i + 1] = (drand48() - 0.5);
        }
    }
    // remove the angular velocity
    zeroDpmsAngularVelocity(dpms);
    // remove the center of mass velocity
    shiftDpmsToVelocity(dpms, 0.0, 0.0);
    // calculate the temperature
    double ke = 0.0;
    int num_vertices = 0;
    for (int id = 0; id < dpms.size(); ++id) {
        dpms[id].calcKineticEnergies();
        ke += dpms[id].kin_eng;
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            num_vertices += 1;
        }
    }
    double temp = 2 * ke / (geomconfig.N_dim * num_vertices * geomconfig.kb);
    // scale the velocities
    double scale_factor = std::sqrt(temp_target / temp);
    std::cout << scale_factor << std::endl;
    std::cout << temp << std::endl;
    std::cout << ke << std::endl;
    for (int id = 0; id < dpms.size(); ++id) {
        for (int i = 0; i < dpms[id].n_vertices; ++i) {
            dpms[id].vel_vertex[2 * i] *= scale_factor;
            dpms[id].vel_vertex[2 * i + 1] *= scale_factor;
        }
    }
}


void compressDpms(std::vector<DPM2D>& dpms, GeomConfig2D& geomconfig, double dr, int N_steps, double phi_target, double damping, int compress_every, std::string dir) {
    // damped compression

    initDataFiles(dir, geomconfig);
    int save_freq = ceil(compress_every / 10);

    // make the logs
    std::ofstream macro_log = createMacroLogFile(dir + "macro_log.csv");
    std::ofstream vertex_log = createVertexLogFile(dir + "vertex_log.csv");
    std::ofstream dpm_log = createDpmLogFile(dir + "dpm_log.csv");
    std::ofstream config_log = createConfigLogFile(dir + "config_log.csv", geomconfig);

    writeMacroConsoleHeader();

    for (int step = 0; step < N_steps; ++step) {
        verletStepDpmList(dpms.size(), dpms, step, damping);

        if (step % save_freq == 0) {
            logDpmList(dpms.size(), dpms, step, save_freq, save_freq * 2, vertex_log, dpm_log, macro_log, config_log, geomconfig);
        }

        if (step % compress_every == 0) {
            writeMacroConsoleHeader();

            double phi = 0.0;
            for (int id = 0; id < dpms.size(); ++id) {
                phi += dpms[id].area;
            }
            phi /= (geomconfig.box_size[0] * geomconfig.box_size[1]);
            if (phi < phi_target) {
                geomconfig.box_size[0] -= dr;
                geomconfig.box_size[1] -= dr;
            }
            else {
                break;
            }
        }
    }

    vertex_log.close();
    dpm_log.close();
    config_log.close();
    macro_log.close();
}


// this is literally the worst garbage code I have ever written - I'm sorry

int getLargestStepNumber(std::string csv_path) {
    std::ifstream file(csv_path);
    std::string line;
    int largest_step = -1;
    if (file.is_open()) {
        while (getline(file, line)) {
            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string column;

            while (getline(iss, column, ',')) {
                columns.push_back(column);
            }

            if (!columns.empty()) {
                try {
                    std::stoi(columns[0]);
                } catch (const std::invalid_argument& ia) {
                    continue;
                }
                int step = std::stoi(columns[0]);
                if (step > largest_step) {
                    largest_step = step;
                }
            }
        }
    }
    file.close();
    return largest_step;
}


void readConfigFileAtStep(std::string config_path, GeomConfig2D& geomconfig, int step) {
    std::ifstream config_file(config_path);
    std::string line;
    if (config_file.is_open()) {
        while (getline(config_file, line)) {
            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string column;

            // Split the line by commas
            while (getline(iss, column, ',')) {
                columns.push_back(column);
            }

            if (!columns.empty()) {
                int step_curr = 0;
                try {
                    step_curr = std::stoi(columns[0]);
                } catch (const std::invalid_argument& ia) {
                    continue;
                }
                if (step_curr == step) {
                    // step,Lx,Ly,dt,kb
                    geomconfig.box_size[0] = std::stod(columns[1]);
                    geomconfig.box_size[1] = std::stod(columns[2]);
                    geomconfig.dt = std::stod(columns[3]);
                    geomconfig.kb = std::stod(columns[4]);
                }
            }
        }
        config_file.close();
    } else {
        std::cerr << "Unable to open config file\n";
    }
}

int getNumberOfDpmsAtStep(std::string dpm_log_path, int step) {
    std::ifstream dpm_log(dpm_log_path);
    std::string line;
    std::vector<int> dpm_ids;
    if (dpm_log.is_open()) {
        while (getline(dpm_log, line)) {
            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string column;

            // Split the line by commas
            while (getline(iss, column, ',')) {
                columns.push_back(column);
            }

            if (!columns.empty()) {
                int step_curr = 0;
                try {
                    step_curr = std::stoi(columns[0]);
                } catch (const std::invalid_argument& ia) {
                    continue;
                }
                if (step_curr == step) {
                    int dpm_id = 0;
                    try {
                        dpm_id = std::stoi(columns[2]);
                    } catch (const std::invalid_argument& ia) {
                        continue;
                    }
                    if (std::find(dpm_ids.begin(), dpm_ids.end(), dpm_id) != dpm_ids.end()) {
                        continue;
                    } else {
                        dpm_ids.push_back(dpm_id);
                    }
                }
            }
        }
        dpm_log.close();
    } else {
        std::cerr << "Unable to open dpm log\n";
    }
    return dpm_ids.size();    
}


void assignVertexDataToDpmAtStep(std::string vertex_log_path, int dpm_id, int step, std::vector<DPM2D>& dpms) {
    int num_vertices = 0;

    // for a given dpm_id, get the vertex data
    std::vector<double> pos_vertex;
    std::vector<double> vel_vertex;
    std::vector<double> force_vertex;
    std::vector<double> acc_vertex;

    std::ifstream vertex_log(vertex_log_path);
    std::string line;
    std::vector<int> vertex_ids;
    if (vertex_log.is_open()) {
        while (getline(vertex_log, line)) {
            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string column;

            // Split the line by commas
            while (getline(iss, column, ',')) {
                columns.push_back(column);
            }

            if (!columns.empty()) {
                int step_curr;
                int dpm_id_curr;
                int vertex_id;
                try {
                    step_curr = std::stoi(columns[0]);
                    dpm_id_curr = std::stoi(columns[2]);
                    vertex_id = std::stoi(columns[3]);
                } catch (const std::invalid_argument& ia) {
                    continue;
                }
                if (step_curr == step && dpm_id_curr == dpm_id) {
                    vertex_ids.push_back(vertex_id);
                    pos_vertex.push_back(std::stod(columns[4]));
                    pos_vertex.push_back(std::stod(columns[5]));
                    vel_vertex.push_back(std::stod(columns[6]));
                    vel_vertex.push_back(std::stod(columns[7]));
                    force_vertex.push_back(std::stod(columns[8]));
                    force_vertex.push_back(std::stod(columns[9]));
                    acc_vertex.push_back(0.0);
                    acc_vertex.push_back(0.0);
                    num_vertices++;
                }
            }
        }
        vertex_log.close();
    } else {
        std::cerr << "Unable to open vertex log\n";
    }

    // assign the stuff to the dpm
    dpms[dpm_id].vel_vertex = vel_vertex;
    dpms[dpm_id].pos_vertex = pos_vertex;
    dpms[dpm_id].force_vertex = force_vertex;
    dpms[dpm_id].acc_vertex = acc_vertex;
    std::cout << "Assigned " << num_vertices << " vertices to dpm " << dpm_id << std::endl;
}


void getForceParamsForDpmAtStep(std::string dpm_log_path, int dpm_id, int step, std::vector<DPM2D>& dpms, ForceParams& forceparams) {
    std::ifstream dpm_log(dpm_log_path);
    std::string line;
    if (dpm_log.is_open()) {
        while (getline(dpm_log, line)) {
            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string column;

            // Split the line by commas
            while (getline(iss, column, ',')) {
                columns.push_back(column);
            }

            if (!columns.empty()) {
                int step_curr;
                int dpm_id;
                try {
                    step_curr = std::stoi(columns[0]);
                    dpm_id = std::stoi(columns[2]);
                } catch (const std::invalid_argument& ia) {
                    continue;
                }
                // check if step = step_curr and dpm_id is not in dpm_ids
                if (step_curr == step && dpm_id) {
                    dpms[dpm_id].pos_dpm[0] = std::stod(columns[3]);
                    dpms[dpm_id].pos_dpm[1] = std::stod(columns[4]);
                    dpms[dpm_id].vel_dpm[0] = std::stod(columns[5]);
                    dpms[dpm_id].vel_dpm[1] = std::stod(columns[6]);
                    dpms[dpm_id].force_dpm[0] = std::stod(columns[7]);
                    dpms[dpm_id].force_dpm[1] = std::stod(columns[8]);
                    dpms[dpm_id].pot_eng = std::stod(columns[9]);
                    dpms[dpm_id].kin_eng = std::stod(columns[10]);
                    dpms[dpm_id].area = std::stod(columns[12]);
                    dpms[dpm_id].perimeter = std::stod(columns[13]);

                    dpms[dpm_id].forceparams.mass_dpm = std::stod(columns[14]);
                    dpms[dpm_id].forceparams.mass_vertex = std::stod(columns[15]);
                    dpms[dpm_id].forceparams.sigma = std::stod(columns[16]);
                    dpms[dpm_id].forceparams.l_0 = std::stod(columns[17]);
                    dpms[dpm_id].forceparams.A_0 = std::stod(columns[18]);
                    dpms[dpm_id].forceparams.theta_0 = std::stod(columns[19]);
                    dpms[dpm_id].forceparams.k = std::stod(columns[20]);
                    dpms[dpm_id].forceparams.k_l = std::stod(columns[21]);
                    dpms[dpm_id].forceparams.k_b = std::stod(columns[22]);
                    dpms[dpm_id].forceparams.k_a = std::stod(columns[23]);
                    dpms[dpm_id].forceparams.damping_coeff = std::stod(columns[24]);
                    dpms[dpm_id].forceparams.Q = std::stod(columns[25]);
                    dpms[dpm_id].forceparams.eta = std::stod(columns[26]);
                }
            }
        }
        dpm_log.close();
    } else {
        std::cerr << "Unable to open vertex log\n";
    }
}

// should probably make it so that we dont have to loop over the entire list every time but oh well - good enough for government work
std::vector<DPM2D> loadDpmData(std::string dir, int step) {
    std::vector<DPM2D> dpms;

    // (if step is -1, then load the last step)
    // get largest step number
    if (step == -1) {
        step = getLargestStepNumber(dir + "dpm_log.csv");
        if (step == -1) {
            std::cerr << "No data found\n";
            return dpms;
        }
    }

    GeomConfig2D geomconfig = GeomConfig2D(0.0, 0.0, 0.0, 0.0);
    ForceParams forceparams = ForceParams(0.0, 0.0, 0.0, 0.0, 0.0);

    // read the config file and define geomconfig
    readConfigFileAtStep(dir + "config_log.csv", geomconfig, step);
    
    // get number of dpms
    int num_dpms = getNumberOfDpmsAtStep(dir + "dpm_log.csv", step);

    // make the dpm list using placeholders
    for (int i = 0; i < num_dpms; i++) {
        dpms.push_back(DPM2D(0.0, 0.0, 0.0, 0, geomconfig, forceparams, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    // for each dpm, assign the vertex data using the vertex log
    for (int i = 0; i < num_dpms; i++) {
        assignVertexDataToDpmAtStep(dir + "vertex_log.csv", i, step, dpms);
        getForceParamsForDpmAtStep(dir + "dpm_log.csv", i, step, dpms, forceparams);
    }
    return dpms;
}