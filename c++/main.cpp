#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "mesh.hpp"
#include "integrals.hpp"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "TrustRegion.hpp"


using namespace std;


int main() {
    //
    //      "|\"
    //      "|_\"
    //
    
    // --------------- load mesh into myMesh -------------------
    // query mesh sizes:
    vector<int> sizes(4, 1);
    std::getMeshLengths(sizes);
    std::cout<<"sizes: ";
    for (int i = 0; i<sizes.size();i++){
      std::cout<<sizes[i]<<" ";
    }
    std::cout<<std::endl;
    // Preallocate arrays in memory: 
    vector<float> Xpoint(sizes[0]);
    vector<float> Ypoint(sizes[0]);
    vector<int> triangles(sizes[1]);
    vector<int> edge1(sizes[2]);
    vector<int> edge2(sizes[3]);
    // Load mesh:
    loadMesh(Xpoint,Ypoint,triangles,edge1,edge2);
    std::mesh myMesh(Xpoint,Ypoint,triangles,edge1,edge2);


    // ----------------------------------------------------


    // Outline:
    // solve for C_0: (K)C = -(f+f_lin)
    // K*C + f + H(c) = 0 (nonlinear system)
  

    // Preallocate arrays: 
    int M = myMesh.getNbNodes();
    int nbBoundaryNodes = myMesh.getNbBoundaryNodes();
    vector<Eigen::Triplet<double>> KmatrixTriplets;
    KmatrixTriplets.reserve(2*9*M);
    vector<Eigen::Triplet<double>> KLinMatrixTriplets;
    KLinMatrixTriplets.reserve(2*9*M+2*3*nbBoundaryNodes);
    Eigen::VectorXd f(2*M);
    Eigen::VectorXd f_lin(2*M);

    //     First integral (K1)
    //     Second integral (K_lin, f_lin)
    //     third integral (K3, f3)
    //  Sparse variables optellen tot eind resultaat
    //      K = K1+K3, f = f3

    integral1(myMesh, KmatrixTriplets); 
    // std::cout<<"KmatrixTriplets: ";
    // //for (int i = 0; i<KmatrixTriplets.size();i++){
    // for (int i = 0; i<10;i++){
    //   std::cout<<KmatrixTriplets[i].value()<<" ";
    // }
    // std::cout<<std::endl;

    integral2lin(myMesh, KLinMatrixTriplets, f_lin);
    // std::cout<<"KLinmatrixTriplets: ";
    // //for (int i = 0; i<KLinMatrixTriplets.size();i++){
    // for (int i = 0; i<10;i++){
    //   std::cout<<KLinMatrixTriplets[i].value()<<" ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"f_lin: ";
    // //for (int i = 0; i<f_lin.size();i++){
    // for (int i = 0; i<10;i++){
    //   std::cout<<f_lin[i]<<" ";
    // }
    // std::cout<<std::endl;

    integral3(myMesh, KmatrixTriplets, f);
    // std::cout<<"KmatrixTriplets: ";
    // //for (int i = 0; i<KmatrixTriplets.size();i++){
    // for (int i = 0; i<10;i++){
    //   std::cout<<KmatrixTriplets[i].value()<<" ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"f: ";
    // //for (int i = 0; i<f.size();i++){
    // for (int i = 0; i<10;i++){
    //   std::cout<<f[i]<<" ";
    // }
    // std::cout<<std::endl;

    // Setup linear system:
    int sizeK = M*2;
    Eigen::SparseMatrix<double> Kmatrix(sizeK,sizeK);
    Kmatrix.setFromTriplets(KmatrixTriplets.begin(), KmatrixTriplets.end());
    Eigen::SparseMatrix<double> KLinMatrix(sizeK,sizeK);
    KLinMatrix.setFromTriplets(KLinMatrixTriplets.begin(), KLinMatrixTriplets.end());

    KLinMatrix += Kmatrix;

    // solve for C: (K)C = -(f+f_lin)
    Eigen::VectorXd newRHS = -(f+f_lin/100);
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.analyzePattern(KLinMatrix);
    solver.factorize(KLinMatrix);
    if(solver.info()!=Eigen::Success) {
        std::cout << "main: error in Eigen Sparse LU factorization KLinMatrix" << std::endl;
    }
    Eigen::VectorXd C0 = solver.solve(newRHS); 
    std::cout<<std::endl<<"Found start solution..."<<std::endl;


    std::ofstream myFile;
    myFile.open("../matlab/cplusplus_output.m");
    myFile<<"C_0 = [ ";
    for (int i = 0; i<C0.size();i++){
      myFile<<C0[i]<<" ";
    }
    myFile<<"];";
    myFile.close();
    std::cout<<"C_0 written to file"<<std::endl;

    //std::cout<<"C0: ";
    //for (int i = 0; i<f.size();i++){
    //for (int i = 0; i<C0.size();i++){
    //  std::cout<<C0[i]<<" ";
    //}
    //std::cout<<std::endl;
    //std::cout<<"End of calculations C0..."<<std::endl;

    // Functie die second integral evalueert voor gegeven C --> H(c)
    Eigen::ArrayXd C0_array = C0.array();
    Eigen::ArrayXd f_array = f.array();
    Eigen::ArrayXd C_array(2*M);
    

    std::trustRegion(myMesh, C_array, std::evaluateCostFunction, C0_array, Kmatrix, f_array);

    std::cout<<std::endl<<"Nonlinear solver converged, writing solution to file..."<<std::endl;

    Eigen::VectorXd C_nonlin = C_array.matrix();

    // write C nonlin to matlab file
    std::ofstream myFileNonLin;
    myFileNonLin.open("../matlab/cnonlin_output.m");
    myFileNonLin<<"C_nonlin = [ ";
    for (int i = 0; i<C_nonlin.size();i++){
      myFileNonLin<<C_nonlin[i]<<" ";
    }
    myFileNonLin<<"];";
    myFileNonLin.close();
    std::cout<<std::endl;

    std::cout<<"File write complete."<<std::endl;

    return 0;
}
