//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/
//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus


#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cassert>
#include "mesh.hpp"

namespace std{

mesh::mesh(std::vector<float> Xpoint, std::vector<float> Ypoint, std::vector<int> triangles, std::vector<int> edge1, std::vector<int> edge2)
    : Xpoint_(Xpoint)
    , Ypoint_(Ypoint)
    , triangles_(triangles)
    // First and last element of edge1_ and edge2_ are the common part!
    , edge1_(edge1)
    , edge2_(edge2)
    {
        assert(Xpoint_.size() == Ypoint_.size());
        assert(triangles_.size() % 3 == 0);
        std::cout<<"Made mesh..."<<std::endl;
    };


const int mesh::getNbElements() {
    return triangles_.size()/3;
}

const int mesh::getNbNodes() {
  return Xpoint_.size();
}

void mesh::getElement(int elementIndex, std::vector<int> &nodeIndices) {
    // return node indices for given element index
    assert (nodeIndices.size() == 3);
    nodeIndices[0] = triangles_[3*elementIndex];
    nodeIndices[1] = triangles_[3*elementIndex+1];
    nodeIndices[2] = triangles_[3*elementIndex+2];
}

void mesh::getNodeCoordinates(int nodeIndex, std::vector<float>& coordinates) {
    assert (coordinates.size() == 2);
    coordinates[0] = Xpoint_[nodeIndex];
    coordinates[1] = Ypoint_[nodeIndex];
}

// Returns number of boundary nodes.
int mesh::getNbBoundaryNodes(bool edge2) {
  if (edge2) {
    return edge2_.size();
  } else {
    return edge1_.size();
  }
}

void mesh::getBoundaryNodes(std::vector<int>& nodeIndices, bool edge2) {
    if (edge2){
        nodeIndices = edge2_;
    }
    else {
        nodeIndices = edge1_;
    }
    return;
}
// end class mesh



// Imports al the numbers from a txt file into a one dimensional std::vector.
template<class T>
void import_matrix_from_txt_file(const char* filename_X, std::vector<T>& v){

    std::ifstream file_X;
    T current;

    file_X.open(filename_X);
    if (file_X.is_open()) {

        // Loop over al the rows of the matrix in the txt file.
        // eof = "end of file."
        int i = 0;
        while (file_X >> current) {
          v[i] = current;
          i++;
        }
        file_X.close();
        std::cout << "Matrix read succes";
    }
    else{
        std::cout << "file open failed";
    }

}


// Load mesh from text files.
void loadMesh(std::vector<float>& Xpoint, std::vector<float>& Ypoint, std::vector<int>& triangles, std::vector<int>& edge1, std::vector<int>& edge2) {
    // loading X-coordinates from txt file
    // std::vector <float> Xpoints;
    import_matrix_from_txt_file("./mesh1/Xpoints.txt", Xpoint);
    std::cout<<"Xpoints: "<<Xpoint[0]<<std::endl;
    // loading Y-coordinates from txt file
    // std::vector <float> Ypoints;
    import_matrix_from_txt_file("./mesh1/Ypoints.txt",Ypoint);
    std::cout<<"Ypoints: "<<Ypoint[0]<<std::endl;
    // loading point indices of triangles from txt file
    // std::vector <int> triangles;
    import_matrix_from_txt_file("./mesh1/triangleLabels.txt",triangles);
    std::cout<<"triangles: "<<triangles[0]<<std::endl;
    // loading point indices of edge from txt file
    // std::vector <int> edge;
    import_matrix_from_txt_file("./mesh1/edge1Labels.txt",edge1);
    std::cout<<"edge1Labels: "<<edge1[0]<<std::endl;
    // loading point indices of edge from txt file
    // std::vector <int> edge;
    import_matrix_from_txt_file("./mesh1/edge2Labels.txt",edge2);
    std::cout<<"edge1Labels: "<<edge2[0]<<std::endl;
}

void getMeshLengths(std::vector<int>& sizes){
  import_matrix_from_txt_file("./mesh1/sizes.txt", sizes);
}



} // end namespace std
