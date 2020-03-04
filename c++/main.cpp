//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/
//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "mesh.hpp"
#include "sibren_functions.hpp"


using namespace std;


int main(){

  // TODO OSCAR: maak hier 1 functie van die alles inlaadt...

    // loading X-coordinates from txt file
    vector <float> Xpoints;
    int xnbpoints=0;
    int xpointRows=0;
    std::import_matrix_from_txt_file("Xpoints.txt",Xpoints,xpointRows,xnbpoints);

    // loading Y-coordinates from txt file
    vector <float> Ypoints;
    int ynbpoints=0;
    int ypointRows=0;
    std::import_matrix_from_txt_file("Ypoints.txt",Ypoints,ypointRows,ynbpoints);
    
    // loading point indices of triangles from txt file
    vector <int> triangles;
    int triangleRows=0;
    int nbtriangles=0;
    std::import_matrix_from_txt_file("triangleLabels.txt",triangles,triangleRows,nbtriangles);

    // loading point indices of edge from txt file
    vector <int> edges;
    int edgeRows=0;
    int nbedge=0;
    std::import_matrix_from_txt_file("edgeLabels.txt",edges,edgeRows,nbedge);

    std::mesh myMesh(Xpoints, Ypoints, triangles, edges);
// --------------------------------------------------------------------

    // (if this code is in your way, put it in the "sem_functions" file)
    // creation of small 1-triangle mesh to test on
    //
    //      |\
    //      |_\
    //
    vector<double> Xpoints{0,0,0.5};
    vector<double> Ypoints{0,1,0};
    vector<int> triangles{1,2,3};
    vector<int> edge{2,3,1}
    mesh testMesh(Xpoints,Ypoints,triangles,edge);
    // -----

  // Outline:
    // K*C + f + H(c) = 0 (nonlinear system)
    // start: solve for C: (K+K_lin)C = -(f+f_lin)

  // Implementation:

    std::vector<Eigen::Triplet<double>> KmatrixTriplets;
    std::vector<Eigen::Triplet<double>> KLinMatrixTriplets;
    KmatrixTriplets.reserve(2*9*myMesh.getNbNodes());
    KLinMatrixTriplets.reserve(2*9*myMesh.getNbNodes()+2*3*myMesh.getNbBoundaryNodes());
    Eigen::vectorXd f(2*myMesh.getNbNodes());
    Eigen::vectorXd f_lin(2*myMesh.getNbNodes());

    //     First integral (K1)
    //     Second integral (K_lin, f_lin)
    //     third integral (K3, f3)
    //  Sparse variables optellen tot eind resultaat
    //      K = K1+K3, f = f3

    integral1(myMesh, KmatrixTriplets);
    integral2lin(myMesh, KLinMatrixTriplets, f_lin);
    integral3(myMesh, KmatrixTriplets, f);

    int sizeK = myMesh.getNbNodes()*2;
    Eigen::SparseMatrix<double,RowMajor> Kmatrix(sizeK,sizeK);
    Kmatrix.setFromTriplets(KmatrixTriplets.begin(), KmatrixTriplets.end());
    Eigen::SparseMatrix<double,RowMajor> KLinMatrix(sizeK,sizeK);
    KLinMatrix.setFromTriplets(KLinMatrixTriplets.begin(), KLinMatrixTriplets.end());

    KLinMatrix += Kmatrix;
    f_lin += f;


    // Solve voor lin oplossing (startwaarde)
      // start: solve for C: (K+K_lin)C = -(f+f_lin)
      Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(KLinMatrix);
      Eigen::VectorXd C0 = chol.sole(f_lin);

    // Functie die second integral evalueert voor gegeven C --> H(c)
    Eigen::VectorXd H(2*myMesh.getNbNodes());
    integral2nonlinear(myMesh, C_current, H);

    // Solve nonlinear createLinearSystem
    // nonlinear_optimise_function(C) = K*C + f + H --> solven voor C tot gelijk aan 0.
    nonlinear_optimise_function(C, K, f, H);


    // Output resultaat C

    // Plot C in Matlab

}
