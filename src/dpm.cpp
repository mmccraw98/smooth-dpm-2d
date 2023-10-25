#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

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

    // update the potential energy
    this->pot_eng_area += this->forceparams.k_a * magnitude_area * magnitude_area / 2.0;
    this->pot_eng += this->pot_eng_area;
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
    this->pot_eng_length += magnitude_length * magnitude_length * this->forceparams.k_l / 2.0;
    this->pot_eng_angle += magnitude_angle * magnitude_angle * this->forceparams.k_b/ 2.0;
    this->pot_eng += this->pot_eng_length + this->pot_eng_angle;
}

void DPM2D::innerDpmForceRoutine() {
    // reset the forces and energies
    this->pot_eng_length = 0.0;
    this->pot_eng_angle = 0.0;
    this->pot_eng_area = 0.0;
    this->pot_eng_int = 0.0;
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
            vertex_log << std::fixed << std::setprecision(precision) << step << "," << this->geomconfig.dt * step << "," << dpm_id << "," << i << "," << this->pos_vertex[this->geomconfig.N_dim * i] << "," << this->pos_vertex[this->geomconfig.N_dim * i + 1] << "," << this->vel_vertex[this->geomconfig.N_dim * i] << "," << this->vel_vertex[this->geomconfig.N_dim * i + 1] << "," << this->force_vertex[this->geomconfig.N_dim * i] << "," << this->force_vertex[this->geomconfig.N_dim * i + 1] << "\n";
        }
    }

    // divide the dpm position by the number of vertices
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) { //
        this->pos_dpm[dim] /= this->n_vertices;
        this->vel_dpm[dim] /= this->n_vertices;
        this->force_dpm[dim] /= this->n_vertices;
    }
    // write the dpm log file
    dpm_log << std::fixed << std::setprecision(precision) << step << "," << this->geomconfig.dt * step << "," << dpm_id << "," << this->pos_dpm[0] << "," << this->pos_dpm[1] << "," << this->vel_dpm[0] << "," << this->vel_dpm[1] << "," << this->force_dpm[0] << "," << this->force_dpm[1] << "," << this->pot_eng_length << "," << this->pot_eng_area << "," << this->pot_eng_angle << "," << this->pot_eng_int << "," << this->pot_eng << "," << this->kin_eng << "," << this->area << "," << this->perimeter << "," << this->forceparams.mass_dpm << "," << this->forceparams.mass_vertex << "," << this->forceparams.sigma << "," << this->forceparams.l_0 << "," << this->forceparams.A_0 << "," << this->forceparams.theta_0 << "," << this->forceparams.k << "," << this->forceparams.k_l << "," << this->forceparams.k_b << "," << this->forceparams.k_a << "\n";
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
    file << "step,t,dpm_id,x,y,vx,vy,fx,fy,pe_length,pe_area,pe_angle,pe_int,pe,ke,area,perim,mass,vertex_mass,vertex_sigma,length_0,area_0,angle_0,k,k_l,k_b,k_a\n";
    return file;
}

void writeMacroConsoleHeader() {
    std::cout << std::string(12 * 7 + 6, '_') << std::endl;
    std::cout << std::setw(12) << "step" << " | " 
              << std::setw(12) << "time" << " | " 
              << std::setw(12) << "pe" << " | " 
              << std::setw(12) << "ke" << " | " 
              << std::setw(12) << "phi" << " | " 
              << std::setw(12) << "temp" << "\n";
    std::cout << std::string(12 * 7 + 6, '_') << std::endl;
}

std::ofstream createMacroLogFile(const std::string& file_path) {
    std::ofstream file;
    file.open(file_path);
    if (!file.is_open()) {
        std::cout << "ERROR: could not open macro log file " << file_path << std::endl;
    }
    file << "step,t,pe,ke,phi,temp\n";

    return file;
}

void writeToMacroLogFile(std::ofstream& macro_log, std::vector<DPM2D>& dpms, int step, GeomConfig2D& geomconfig, int precision, int console_log_freq) {
    if (!macro_log) {
        std::cout << "ERROR: could not write to macro log file" << std::endl;
    }
    // calculate the macro variables
    double pe = 0.0;
    double ke = 0.0;
    double phi = 0.0;
    for (int id = 0; id < dpms.size(); ++id) {
        pe += dpms[id].pot_eng;
        ke += dpms[id].kin_eng;
        phi += dpms[id].area;
    }
    // calculate the temperature
    double temp = 2 * ke / (geomconfig.N_dim * dpms.size());
    // calculate the packing fraction
    phi /= (geomconfig.box_size[0] * geomconfig.box_size[1]);

    // write the macro log file
    macro_log << std::fixed << std::setprecision(precision) << step << "," << dpms[0].geomconfig.dt * step << "," << pe << "," << ke << "," << phi << "," << temp << "\n";
    // log to console the same thing
    if (step % console_log_freq == 0) {
       std::cout << std::setw(12) << std::fixed << std::setprecision(3) << step << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << (dpms[0].geomconfig.dt * step) << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << pe << " | " 
          << std::setw(12) << std::scientific << std::setprecision(3) << ke << " | " 
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
    std::string command = "mkdir " + dir;
    system(command.c_str());

    // create the config file in the dir
    std::ofstream config_file;
    config_file.open(dir + "config.csv");
    config_file << "Lx,Ly,dt,kb" << std::endl;
    config_file << geomconfig.box_size[0] << "," << geomconfig.box_size[1] << "," << geomconfig.dt << "," << geomconfig.dt << std::endl;
    config_file.close();
}

void logDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step, int save_freq, int console_log_freq, std::ofstream& vertex_log, std::ofstream& dpm_log, std::ofstream& macro_log, GeomConfig2D& geomconfig) {
    if (step % save_freq == 0) {
        for (int id = 0; id < num_dpms; ++id) {
            dpms[id].writeToLogFiles(vertex_log, dpm_log, step, id);
        }

        // write the macro variables to the log file
        writeToMacroLogFile(macro_log, dpms, step, geomconfig, 10, console_log_freq);
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
                other_dpm.pot_eng_int += this->forceparams.k * magnitude * magnitude / 2.0;
                this->pot_eng_int += this->forceparams.k * magnitude * magnitude / 2.0;
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
                other_dpm.pot_eng_int += this->forceparams.k * magnitude * magnitude / 2.0;
                this->pot_eng_int += this->forceparams.k * magnitude * magnitude / 2.0;
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
        this->setParticleSegmentForceEnergy(other_dpm, other_vertex_i);
        // calculate the particle - vertex interaction
        this->setParticleVertexForceEnergy(other_dpm, other_vertex_i);
    }
}

void verletStepDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step) {
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

    for (int id = 0; id < num_dpms; ++id) {  // TODO parallelize this
        dpms[id].verletVelocityStep();
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


double setLinearHarmonicForces(std::vector<double>& force, std::vector<double>& pos, std::vector<double>& distance_vector, int num_disks, GeomConfig2D& geomconfig, ForceParams& forceparams) {
    // calculate the area of all the particles
    double area = 0.0;
    
    // define a force routine
    for (int i = 0; i < num_disks; ++i) {
        area += M_PI * forceparams.sigma * forceparams.sigma / 4.0;
        for (int j = i + 1; j < num_disks; ++j) {
            // calculate the distance between the two disks
            for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
                distance_vector[dim] = getPbcDistPoint(pos[i * geomconfig.N_dim + dim], pos[j * geomconfig.N_dim + dim], geomconfig.box_size[dim]);
            }
            double length = getVectLength(distance_vector, geomconfig);
            double overlap = forceparams.sigma - length;  // NOTE - in bidisperse case, this will be (sigma_i + sigma_j) / 2 - length
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


double adamMinimizeDiskForces(std::vector<double>& pos, GeomConfig2D& geomconfig, ForceParams& forceparams, int num_disks, int max_steps, double alpha, double beta1, double beta2, double epsilon) {
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
        area = setLinearHarmonicForces(gradient, pos, distance_vector, num_disks, geomconfig, forceparams);
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