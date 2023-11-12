#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

#include "misc.hpp"

void writeMacroConsoleHeader() {  // TODO make this more general
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

std::vector<double> getCircleCoords(double cx, double cy, double radius, double n_vertices) {
    // reserve a vector of size 2*n_vertices and fill it with the x and y coordinates of the vertices
    std::vector<double> coords(2*n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        coords[2 * i] = cx + radius * cos(2 * M_PI * i / n_vertices);
        coords[2 * i + 1] = cy + radius * sin(2 * M_PI * i / n_vertices);
    }
    return coords;
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

std::vector<double> generateRandomUniformVector(int N, double min, double max, double seed) {
    std::vector<double> randomUniformVector(N, 0.0);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(min, max);
    for (int i = 0; i < N; ++i) {
        randomUniformVector[i] = dist(gen);
    }
    return randomUniformVector;
}

std::vector<double> generateRandomNormalVector(int N, double mean, double std_deviation, double seed) {
    std::vector<double> randomNormalVector(N, 0.0);
    std::mt19937 gen(seed);
    std::normal_distribution<double> dist(mean, std_deviation);
    for (int i = 0; i < N; ++i) {
        randomNormalVector[i] = dist(gen);
    }
    return randomNormalVector;
}