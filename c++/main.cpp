#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "mesh.hpp"
#include "integrals.hpp"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "LM.hpp"


using namespace std;


int main() {
    // (if this code is in your way, put it in the "sem_functions" file)
    // creation of small 1-triangle mesh to test on
    //
    //      "|\"
    //      "|_\"
    //
    //vector<float> Xpoint{0,0,0.5};
    //vector<float> Ypoint{0,1,0};
    //vector<int> triangles{1,2,3};
    //vector<int> edge{2,3,1};
    //std::mesh testMesh(Xpoint,Ypoint,triangles,edge);
    // -----

    // load mesh into variable   myMesh -------------------
    // vector<float> Xpoints, Ypoints;
    // vector<int> triangles, edge;
    vector<int> sizes(4, 1);
    // sizes = [points.length, triangles.length, edge1.length, edge2.length]
    std::getMeshLengths(sizes);
    std::cout<<"sizes: ";
    for (int i = 0; i<sizes.size();i++){
      std::cout<<sizes[i]<<" ";
    }
    std::cout<<std::endl;
    vector<float> Xpoint(sizes[0]);
    vector<float> Ypoint(sizes[0]);
    vector<int> triangles(sizes[1]);
    vector<int> edge1(sizes[2]);
    vector<int> edge2(sizes[3]);
    loadMesh(Xpoint,Ypoint,triangles,edge1,edge2);

    std::mesh myMesh(Xpoint,Ypoint,triangles,edge1,edge2);
    // ----------------------------------------------------


    // Outline:
    // K*C + f + H(c) = 0 (nonlinear system)
    // start: solve for C: (K+K_lin)C = -(f+f_lin)

    // Implementation:
    int M = myMesh.getNbNodes();
    int nbBoundaryNodes = myMesh.getNbBoundaryNodes();
    std::cout<<"been here"<<std::endl;
    vector<Eigen::Triplet<double>> KmatrixTriplets;
    KmatrixTriplets.reserve(2*9*M);
    vector<Eigen::Triplet<double>> KLinMatrixTriplets;
    KLinMatrixTriplets.reserve(2*9*M+2*3*nbBoundaryNodes);
    //KmatrixTriplets.reserve(2*9*M);
    //KLinMatrixTriplets.reserve(2*9*M+2*3*nbBoundaryNodes);
    Eigen::VectorXd f(2*M);
    Eigen::VectorXd f_lin(2*M);

    //     First integral (K1)
    //     Second integral (K_lin, f_lin)
    //     third integral (K3, f3)
    //  Sparse variables optellen tot eind resultaat
    //      K = K1+K3, f = f3
    integral1(myMesh, KmatrixTriplets); // Sibren_functions OK!
    std::cout<<"KmatrixTriplets: ";
    //for (int i = 0; i<KmatrixTriplets.size();i++){
    for (int i = 0; i<10;i++){
      std::cout<<KmatrixTriplets[i].value()<<" ";
    }
    std::cout<<std::endl;
    integral2lin(myMesh, KLinMatrixTriplets, f_lin);
    std::cout<<"KLinmatrixTriplets: ";
    //for (int i = 0; i<KLinMatrixTriplets.size();i++){
    for (int i = 0; i<10;i++){
      std::cout<<KLinMatrixTriplets[i].value()<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"f_lin: ";
    //for (int i = 0; i<f_lin.size();i++){
    for (int i = 0; i<10;i++){
      std::cout<<f_lin[i]<<" ";
    }
    std::cout<<std::endl;
    integral3(myMesh, KmatrixTriplets, f);
    std::cout<<"KmatrixTriplets: ";
    //for (int i = 0; i<KmatrixTriplets.size();i++){
    for (int i = 0; i<10;i++){
      std::cout<<KmatrixTriplets[i].value()<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"f: ";
    //for (int i = 0; i<f.size();i++){
    for (int i = 0; i<10;i++){
      std::cout<<f[i]<<" ";
    }
    std::cout<<std::endl;


    int sizeK = M*2;
    Eigen::SparseMatrix<double> Kmatrix(sizeK,sizeK);
    Kmatrix.setFromTriplets(KmatrixTriplets.begin(), KmatrixTriplets.end());
    Eigen::SparseMatrix<double> KLinMatrix(sizeK,sizeK);
    KLinMatrix.setFromTriplets(KLinMatrixTriplets.begin(), KLinMatrixTriplets.end());

    KLinMatrix += Kmatrix;


    // Solve voor lin oplossing (startwaarde)
    // start: solve for C: (K+K_lin)C = -(f+f_lin)

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(KLinMatrix);
    Eigen::VectorXd C0 = chol.solve(-f_lin);

    std::ofstream myFile;
    myFile.open("../matlab/cplusplus_output.m");
    myFile<<"C_0 = [ ";
    for (int i = 0; i<C0.size();i++){
      myFile<<C0[i]<<" ";
    }
    myFile<<"]';";
    std::cout<<std::endl;
    std::cout<<"C0: ";
    //for (int i = 0; i<f.size();i++){
    for (int i = 0; i<10;i++){
      std::cout<<f[i]<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"End of calculations C0..."<<std::endl;

    // Functie die second integral evalueert voor gegeven C --> H(c)

    Eigen::ArrayXd C0_array = C0.array();
    bool flag = (C0_array.maxCoeff() < 1e6 && C0_array.minCoeff() > -1e6);
    std::cout<<"flag: "<<flag;
    
    Eigen::ArrayXd C_array(2*M);

    std::minimize_lm(myMesh, C_array, std::integral2nonlinear, C0_array);



    // integral2nonlinear(myMesh, C_current, H); // Sibren_functions OK!


    // Solve nonlinear createLinearSystem
    // nonlinear_optimise_function(C) = K*C + f + H --> solven voor C tot gelijk aan 0.
    // nonlinear_optimise_function(C, K, f, H);


    // Output resultaat C

    // Plot C in Matlab

    return 0;
}
