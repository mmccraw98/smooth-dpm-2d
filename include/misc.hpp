#ifndef MISC_HPP
#define MISC_HPP

#include <vector>

std::vector<double> getCircleCoords(double cx, double cy, double radius, double n_vertices);
std::vector<double> generateLatticeCoordinates(int N, double lx, double ly);
std::vector<double> generateRandomUniformVector(int N, double min, double max, double seed);
std::vector<double> generateRandomNormalVector(int N, double mean, double std_deviation, double seed);

#endif // MISC_HPP