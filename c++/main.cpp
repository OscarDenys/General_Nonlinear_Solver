#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "mesh.hpp"
#include "integrals.hpp"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"


using namespace std;


int main() {
    // (if this code is in your way, put it in the "sem_functions" file)
    // creation of small 1-triangle mesh to test on
    //
    //      "|\"
    //      "|_\"
    //
    vector<float> Xpoint{0,0,0.5};
    vector<float> Ypoint{0,1,0};
    vector<int> triangles{1,2,3};
    vector<int> edge{2,3,1};
    std::mesh testMesh(Xpoint,Ypoint,triangles,edge);
    // -----

    // load mesh into variable   myMesh -------------------
    // vector<float> Xpoints, Ypoints;
    // vector<int> triangles, edge;
    loadMesh(Xpoint,Ypoint,triangles,edge);
    std::mesh myMesh(Xpoint,Ypoint,triangles,edge);
    // ----------------------------------------------------


    // Outline:
    // K*C + f + H(c) = 0 (nonlinear system)
    // start: solve for C: (K+K_lin)C = -(f+f_lin)

    // Implementation:

    vector<Eigen::Triplet<double>> KmatrixTriplets;
    vector<Eigen::Triplet<double>> KLinMatrixTriplets;
    int M = myMesh.getNbNodes();
    int nbBoundaryNodes = myMesh.getNbBoundaryNodes();
    KmatrixTriplets.reserve(2*9*M);
    KLinMatrixTriplets.reserve(2*9*M+2*3*nbBoundaryNodes);
    Eigen::VectorXd f(2*M);
    Eigen::VectorXd f_lin(2*M);

    //     First integral (K1)
    //     Second integral (K_lin, f_lin)
    //     third integral (K3, f3)
    //  Sparse variables optellen tot eind resultaat
    //      K = K1+K3, f = f3

    integral1(myMesh, KmatrixTriplets); // Sibren_functions OK!
    integral2lin(myMesh, KLinMatrixTriplets, f_lin);
    integral3(myMesh, KmatrixTriplets, f);

    int sizeK = M*2;
    Eigen::SparseMatrix<double> Kmatrix(sizeK,sizeK);
    Kmatrix.setFromTriplets(KmatrixTriplets.begin(), KmatrixTriplets.end());
    Eigen::SparseMatrix<double> KLinMatrix(sizeK,sizeK);
    KLinMatrix.setFromTriplets(KLinMatrixTriplets.begin(), KLinMatrixTriplets.end());

    KLinMatrix += Kmatrix;
    f_lin += f;


    // Solve voor lin oplossing (startwaarde)
    // start: solve for C: (K+K_lin)C = -(f+f_lin)

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(KLinMatrix);
    Eigen::VectorXd C0 = chol.solve(f_lin);

    // Functie die second integral evalueert voor gegeven C --> H(c)
    Eigen::VectorXd H(2*M);
    // integral2nonlinear(myMesh, C_current, H); // Sibren_functions OK!


    // Solve nonlinear createLinearSystem
    // nonlinear_optimise_function(C) = K*C + f + H --> solven voor C tot gelijk aan 0.
    // nonlinear_optimise_function(C, K, f, H);


    // Output resultaat C

    // Plot C in Matlab

    return 0;
}
