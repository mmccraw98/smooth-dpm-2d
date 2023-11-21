#ifndef MISC_HPP
#define MISC_HPP

#include <vector>

std::vector<double> getCircleCoords(double cx, double cy, double radius, double n_vertices);
std::vector<double> generateLatticeCoordinates(int N, double lx, double ly);
std::vector<double> generateRandomUniformVector(int N, double min, double max, double seed);
std::vector<double> generateRandomNormalVector(int N, double mean, double std_deviation, double seed);

inline double generateRandomNormal() {
    // Static to initialize only once
    static std::random_device rd;
    static std::default_random_engine generator(rd());
    static std::normal_distribution<double> distribution(0.0, 1.0);

    return distribution(generator);
}

#endif // MISC_HPP