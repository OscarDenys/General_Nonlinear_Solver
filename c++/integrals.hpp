#ifndef integrals_hpp
#define integrals_hpp
#include "mesh.hpp"
#include "Eigen/SparseCore"
#include <vector>

typedef Eigen::Triplet<double> Trip;

namespace std {

    void integral1(mesh const &myMesh, std::vector<Trip> &K);

    void integral2lin(mesh &myMesh, vector<Trip> & K_lin, Eigen::VectorXd & f_lin);
    void integral2nonlinear(mesh & myMesh, Eigen::VectorXd &H, Eigen::VectorXd &C);

    void integral3(mesh &myMesh, vector<Trip> & K, Eigen::VectorXd & f);

} // namespace std

#endif