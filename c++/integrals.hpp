#ifndef integrals_hpp
#define integrals_hpp
#include "mesh.hpp"
#include "Eigen/SparseCore"
#include <vector>

typedef Eigen::Triplet<double> Trip;

namespace std {

    void integral1(std::mesh &myMesh, std::vector<Trip> &K);

    void integral2lin(mesh &myMesh, vector<Trip> & K_lin, Eigen::VectorXd & f_lin);
    Eigen::ArrayXd integral2nonlinear(Eigen::ArrayXd &C, std::mesh & myMesh);

    void integral3(mesh &myMesh, vector<Trip> & K, Eigen::VectorXd & f);

    double detJac(vector<float> P1, vector<float> P2, vector<float> P3);
    void applySigmaAndAddCommonPart(std::vector<double> &result, std::vector<double> const &commonPart);

    void evaluateRespiration(int nodeIndex, Eigen::ArrayXd &prevSol, double &Ru, double &Rv);
    void evaluateRespiration(int nodeIndex1, int nodeIndex2, Eigen::ArrayXd &prevSol, double &Ru, double &Rv);
    double evaluateRu(double Cu, double Cv);
    double evaluateRv(double Cu, double Cv, double Ru);
    void evaluateCostFunction( Eigen::SparseMatrix<double> &K, Eigen::ArrayXd &f, Eigen::ArrayXd &C, Eigen::ArrayXd &F, std::mesh &myMesh);

} // namespace std

#endif
