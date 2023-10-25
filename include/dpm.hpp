#ifndef DPM_HPP
#define DPM_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

struct GeomConfig2D {
    int box_size[2];  // Store x_length and y_length
    int N_dim;  // number of dimensions
    double dt;  // time step
    double kb;  // boltzmann constant

    GeomConfig2D(double x_length, double y_length, double dt, double kb) {
        box_size[0] = x_length;
        box_size[1] = y_length;
        N_dim = 2;
        this->dt = dt;
        this->kb = kb;
    }
};

struct ForceParams {
    double k;  // particle - particle interaction
    double k_l;  // bond stretching
    double k_b;  // bond bending
    double k_a;  // area constraint
    double mass_vertex;  // mass of the vertices
    double mass_dpm;  // mass of the dpm
    double sigma;  // particle - particle interaction distance (width of particles)
    double l_0;  // bond length
    double theta_0;  // bond angle
    double A_0;  // preferred area
    double damping_coeff;  // damping coefficient

    ForceParams(double k, double k_l, double k_b, double k_a, double mass_vertex) {
        this->k = k;
        this->k_l = k_l;
        this->k_b = k_b;
        this->k_a = k_a;
        this->mass_vertex = mass_vertex;
        this->damping_coeff = 0.0;  // default to be 0!
        // these are default values and should be overwritten when the DPM2D is initiated
        double sigma = 1.0;  // particle - particle interaction distance (width of particles)
        double l_0 = 1.0;  // bond length
        double theta_0 = M_PI;  // bond angle
        double A_0 = 1.0;  // preferred area
        double mass_dpm = 1.0;  // mass of the dpm
    }
};

class DPM2D {
    public:
        // Member variables
        std::vector<double> pos_vertex, vel_vertex, force_vertex, acc_vertex;
        std::vector<double> bond_lengths, bond_angles;
        std::vector<double> pos_dpm, vel_dpm, force_dpm;
        std::vector<double> pos_i, pos_j, pos_k;
        std::vector<double> l_ij, l_ki;
        std::vector<bool> vertex_is_active;
        double pot_eng, kin_eng, pot_eng_length, pot_eng_angle, pot_eng_area, pot_eng_int;
        int n_vertices;
        double area, perimeter;
        ForceParams forceparams;
        GeomConfig2D& geomconfig;

        // Member functions
        void initializeVectors();
        inline void verletPositionStep();
        inline void verletVelocityStep();
        inline void setDpmPosition();
        inline void assignRandomNormalVelocities(const double& sigma, const double& seed);
        inline void centerVelocities(double vx, double vy);
        inline void setDpmVelocity();
        inline void calcKineticEnergies();
        inline double getTemperature();
        inline double rescaleVelocitiesToTemp(double T_target);
        inline std::vector<double> getVertexPos(const int i);
        inline void setPos_i(const int i);
        inline void setPos_j(const int j);
        inline void setPos_k(const int k);
        inline int getNextVertex(const int i);
        inline int getPrevVertex(const int i);
        inline void setProjLengthToPos_k(const double proj_length, const double bond_length);
        void setAreaForceEnergy(const int i, const int k);
        void setBondBendStretchForcesEnergies(const int i, const int j, const int k);
        void innerDpmForceRoutine();
        void calcBondLengthsAnglesAreaPerimeter();
        void writeToLogFiles(std::ofstream& vertex_log, std::ofstream& dpm_log, int step, int dpm_id, int precision=10);
        void setParticleVertexForceEnergy(DPM2D& other_dpm, int other_vertex_i);
        void setParticleSegmentForceEnergy(DPM2D& other_dpm, int other_vertex_i);
        void setInteractionForceEnergy(DPM2D& other_dpm);
        void innerDpmForceRoutinePlateCompression(std::vector<double> wall_bounds, double wall_strength, std::vector<double>& force_area, std::vector<double>& boundary_pos);
        void adamMinimizeDpmPlateForces(int max_steps, double alpha, double beta1, double beta2, double epsilon, std::vector<double> wall_bounds, double wall_strength);
        void gradDescMinDpmPlateForces(int max_steps, double eta, double tol, std::vector<double> wall_bounds, double wall_strength, std::vector<double>& force_area, std::vector<double>& boundary_pos);

        // Default Constructor
        DPM2D(double cx, double cy, double radius, double n_vertices, GeomConfig2D& geomconfig, ForceParams forceparams, double length_diam_ratio = 2.0, double vx=0.0, double vy=0.0, double T_0=0.0, double seed=12345.0);
        
        // Default Destructor
        ~DPM2D();
};

inline void DPM2D::verletPositionStep() {
    // update the positions
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            this->pos_vertex[this->geomconfig.N_dim * i + dim] += this->vel_vertex[this->geomconfig.N_dim * i + dim] * this->geomconfig.dt + 0.5 * this->force_vertex[this->geomconfig.N_dim * i + dim] * this->geomconfig.dt * this->geomconfig.dt / this->forceparams.mass_vertex;
        }
    }
}

inline void DPM2D::verletVelocityStep() {
    // update the velocities and accelerations
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            this->vel_vertex[this->geomconfig.N_dim * i + dim] += 0.5 * this->geomconfig.dt * (this->force_vertex[this->geomconfig.N_dim * i + dim] / this->forceparams.mass_vertex + this->acc_vertex[this->geomconfig.N_dim * i + dim]);
            this->acc_vertex[this->geomconfig.N_dim * i + dim] = this->force_vertex[this->geomconfig.N_dim * i + dim] / this->forceparams.mass_vertex;
        }
    }
}

inline void DPM2D::setDpmPosition() {
    // sum over vertices and calculate the average position
    std::fill(this->pos_dpm.begin(), this->pos_dpm.end(), 0.0);
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            this->pos_dpm[dim] += this->pos_vertex[this->geomconfig.N_dim * i + dim];
        }
    }
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->pos_dpm[dim] /= this->n_vertices;
    }
}

inline void DPM2D::assignRandomNormalVelocities(const double& sigma, const double& seed) {
    // pick the velocities from random numbers that depend on an input seed variable
    // get number of points in velocities array
    std::fill(this->vel_vertex.begin(), this->vel_vertex.end(), 0.0);
    srand48(seed);
    for (int i = 0; i < this->n_vertices; ++i) {
        this->vel_vertex[this->geomconfig.N_dim * i] = sigma * (drand48() - 0.5);
        this->vel_vertex[this->geomconfig.N_dim * i + 1] = sigma * (drand48() - 0.5);
    }
}

inline void DPM2D::setDpmVelocity() {
    // sum over vertices and calculate the average velocity
    std::fill(this->vel_dpm.begin(), this->vel_dpm.end(), 0.0);
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            this->vel_dpm[dim] += this->vel_vertex[this->geomconfig.N_dim * i + dim];
        }
    }
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->vel_dpm[dim] /= this->n_vertices;
    }
}

inline void DPM2D::centerVelocities(double vx, double vy) {
    this->setDpmVelocity();
    // subtract the average velocity from each velocity
    for (int i = 0; i < this->n_vertices; ++i) {
        this->vel_vertex[this->geomconfig.N_dim * i] -= this->vel_dpm[0] - vx;
        this->vel_vertex[this->geomconfig.N_dim * i + 1] -= this->vel_dpm[1] - vy;
    }
    // recalculate the dpm velocity
    this->setDpmVelocity();
}

inline void DPM2D::calcKineticEnergies() {
    // calculate the kinetic energy of the vertices
    this->kin_eng = 0.0;
    for (int i = 0; i < this->n_vertices; ++i) {
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            this->kin_eng += this->forceparams.mass_vertex * this->vel_vertex[this->geomconfig.N_dim * i + dim] * this->vel_vertex[this->geomconfig.N_dim * i + dim] / 2.0;
        }
    }
    // calculate the kinetic energy of the dpm
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->kin_eng += this->forceparams.mass_dpm * this->vel_dpm[dim] * this->vel_dpm[dim] / 2.0;
    }
}

inline double DPM2D::getTemperature() {
    return 2.0 * this->kin_eng / (this->n_vertices * this->forceparams.mass_vertex * this->geomconfig.kb);
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
        for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
            this->vel_vertex[this->geomconfig.N_dim * i + dim] *= scale_factor;
        }
    }
    this->setDpmVelocity();
    this->calcKineticEnergies();
    return this->getTemperature();
}

inline void DPM2D::setPos_i(const int i) {
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->pos_i[dim] = this->pos_vertex[this->geomconfig.N_dim * i + dim];
    }
}

inline void DPM2D::setPos_j(const int j) {
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->pos_j[dim] = this->pos_vertex[this->geomconfig.N_dim * j + dim];
    }
}

inline void DPM2D::setPos_k(const int k) {
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        this->pos_k[dim] = this->pos_vertex[this->geomconfig.N_dim * k + dim];
    }
}

inline void DPM2D::setProjLengthToPos_k(const double proj_length, const double bond_length) {  // requires that pos_j is set
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
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
    std::vector<double> pos_i(this->geomconfig.N_dim);
    for (int dim = 0; dim < this->geomconfig.N_dim; ++dim) {
        pos_i[dim] = this->pos_vertex[this->geomconfig.N_dim * i + dim];
    }
    return pos_i;
}

inline double getDotProd(const std::vector<double>& vect1, const std::vector<double>& vect2, const GeomConfig2D& geomconfig) {
    double dot_prod = 0.0;
    for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
        dot_prod += vect1[dim] * vect2[dim];
    };
    return dot_prod;
}

inline double getPbcDistPoint(const double& x1, const double& x2, const double& axis_length) {
    double dx = x1 - x2;
    dx -= axis_length * std::round(dx / axis_length);
    return dx;
}

inline double getDist(const double *this_pos, const double *other_pos, const GeomConfig2D& geomconfig) {
    double dist_sq = 0.0;
    double delta = 0.0;
    for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
        delta = getPbcDistPoint(this_pos[dim], other_pos[dim], geomconfig.box_size[dim]);
        dist_sq += delta * delta;
    };
    return std::sqrt(dist_sq);
}

inline double getVectLength(const std::vector<double>& vect, const GeomConfig2D& geomconfig) {
    double vect_length_sq = 0.0;
    for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
        vect_length_sq += vect[dim] * vect[dim];
    }
    return std::sqrt(vect_length_sq);
}

inline void setDistVect(std::vector<double>& dist_vect, const std::vector<double>& this_pos, const std::vector<double>& other_pos, const GeomConfig2D& geomconfig) {
    for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
        dist_vect[dim] = getPbcDistPoint(this_pos[dim], other_pos[dim], static_cast<double>(geomconfig.box_size[dim]));
    }
}

inline std::vector<double> getDistVect(const std::vector<double>& this_pos, const std::vector<double>& other_pos, const GeomConfig2D& geomconfig) {
    std::vector<double> dist_vect(geomconfig.N_dim);
    for (int dim = 0; dim < geomconfig.N_dim; ++dim) {
        dist_vect[dim] = getPbcDistPoint(this_pos[dim], other_pos[dim], static_cast<double>(geomconfig.box_size[dim]));
    }
    return dist_vect;
}

inline double calcProjection(const std::vector<double> this_pos, const std::vector<double> segment_origin_pos, const std::vector<double> segment_vec, const double& segment_length, const GeomConfig2D& geomconfig) {
    return (getPbcDistPoint(this_pos[0], segment_origin_pos[0], geomconfig.box_size[0]) * segment_vec[0] + getPbcDistPoint(this_pos[1], segment_origin_pos[1], geomconfig.box_size[1]) * segment_vec[1]) / segment_length;
}

std::vector<double> getCircleCoords(double cx, double cy, double radius, double n_vertices);
std::ofstream createDpmLogFile(const std::string& file_path);
std::ofstream createVertexLogFile(const std::string& file_path);
std::ofstream createMacroLogFile(const std::string& file_path);
void initDataFiles(std::string dir, GeomConfig2D geomconfig);
void writeToMacroLogFile(std::ofstream& macro_log, std::vector<DPM2D>& dpms, int step, GeomConfig2D& geomconfig, int precision, int console_log_freq);
void writeMacroConsoleHeader();
void tracerTest(std::string dir, int num_vertices);
std::vector<double> generateLatticeCoordinates(int N, double lx, double ly);
void verletStepDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step);
void logDpmList(int num_dpms, std::vector<DPM2D>& dpms, int step, int save_freq, int console_log_freq, std::ofstream& vertex_log, std::ofstream& dpm_log, std::ofstream& macro_log, GeomConfig2D& geomconfig);
double setLinearHarmonicForces(std::vector<double>& force, std::vector<double>& pos, std::vector<double>& distance_vector, int num_disks, GeomConfig2D& geomconfig, ForceParams& forceparams);
double adamMinimizeDiskForces(std::vector<double>& pos, GeomConfig2D& geomconfig, ForceParams& forceparams, int num_disks, int max_steps, double alpha, double beta1, double beta2, double epsilon);
#endif // DPM_HPP