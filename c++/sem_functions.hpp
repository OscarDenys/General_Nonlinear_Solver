#ifndef sem_functions_hpp
#define sem_functions_hpp
#include "constants.hpp"
#include "mesh.hpp"
#include <Eigen/SparseCore>
#include <cassert>
#include <cmath>
#include <vector>

typedef Eigen::Triplet<double> Trip;

namespace std {

    double detJac(vector<float> P1, vector<float> P2, vector<float> P3);

    void evaluateRespiration(int nodeIndex, std::vector<double> prevSol, double Ru, double Rv);
    void evaluateRespiration(int nodeIndex1, int nodeIndex2, std::vector<double> prevSol, double Ru, double Rv);
    double evaluateRu(double Cu, double Cv);
    double evaluateRv(double Cu, double Cv, double Ru);

    void integral2lin(mesh & myMesh, vector<Trip> & K_lin, Eigen::VectorXd & f_lin);

    void integral3(mesh & myMesh, vector<Trip> & K, Eigen::VectorXd & f);
}


#endif	