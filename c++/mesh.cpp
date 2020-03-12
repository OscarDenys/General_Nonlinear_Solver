//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/
//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus


#include "mesh.hpp"
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cassert>

namespace meshClass{


class mesh{

    private:
        const std::vector<float> Xpoint_;
        const std::vector<float> Ypoint_;

        // Structure triangles
        //  - one dimensional std::vector containing the indices of the Xpoints and Ypoints std::vectors.
        //  - each pair of three elements (index%3 = {0,1,2}) forms a triangle.
        std::vector<int> triangles_;

        // Structure edge
        //  - one dimensional vector containing the indices of the Xpoint and Ypoint vectors.
        //  - the indices are given in clockwise direction around the edge (starting from the top of the pear).

        std::vector<int> edge1_;
        std::vector<int> edge2_;


    public:
    mesh(std::vector<float> Xpoint, std::vector<float> Ypoint, std::vector<int> triangles, std::vector<int> edge1, std::vector<int> edge2)
    : Xpoint_(Xpoint)
    , Ypoint_(Ypoint)
    , triangles_(triangles)
    // First and last element of edge1_ and edge2_ are the common part!
    , edge1_(edge1)
    , edge2_(edge2)
    {
        assert(Xpoint_.size() == Ypoint_.size());
        assert(triangles_.size() % 3 == 0);
    }


    const int getNbElements() {
        return triangles_.size()/3;
    }

    const int getNbNodes() {
      return Xpoint_.size();
    }

    const void getElement(int elementIndex, std::vector<int> &nodeIndices) {
        // return node indices for given element index
        assert (nodeIndices.size() == 3);
        nodeIndices[0] = triangles_[3*elementIndex];
        nodeIndices[1] = triangles_[3*elementIndex+1];
        nodeIndices[2] = triangles_[3*elementIndex+2];
    }

    const void getNodeCoordinates(int nodeIndex, std::vector<float> coordinates) {
        assert (coordinates.size() == 2);
        coordinates[0] = Xpoint_[nodeIndex];
        coordinates[1] = Ypoint_[nodeIndex];
    }

    // Returns number of boundary nodes.
    const int getNbBoundaryNodes(bool edge2 = true) {
      if (edge2) {
        return edge2_.size();
      } else {
        return edge1_.size();
      }
    }

    const void getBoundaryNodes(std::vector<int>& nodeIndices, bool edge2 = true) {
        if (edge2){
            nodeIndices = edge2_;
        }
        else {
            nodeIndices = edge1_;
        }
        return;
    }




}; // end class mesh



// Imports al the numbers from a txt file into a one dimensional std::vector.
template<class T>
void import_matrix_from_txt_file(const char* filename_X, std::vector <T>& v){

    std::ifstream file_X;
    int current;

    file_X.open(filename_X);
    if (file_X.is_open()) {

        // Loop over al the rows of the matrix in the txt file.
        // eof = "end of file."
        while (file_X >> current) {
          v.push_back(current);
        }
        file_X.close();
        std::cout << "Matrix read succes";
    }
    else{
        std::cout << "file open failed";
    }

}


// Load mesh from text files.
void loadMesh(std::vector <float> Xpoint, std::vector <float> Ypoint, std::vector <int> triangles, std::vector <int> edge1, std::vector<int> edge2) {
    // loading X-coordinates from txt file
    // std::vector <float> Xpoints;
    meshClass::import_matrix_from_txt_file("Xpoints.txt", Xpoint);

    // loading Y-coordinates from txt file
    // std::vector <float> Ypoints;
    meshClass::import_matrix_from_txt_file("Ypoints.txt",Ypoint);

    // loading point indices of triangles from txt file
    // std::vector <int> triangles;
    meshClass::import_matrix_from_txt_file("triangleLabels.txt",triangles);

    // loading point indices of edge from txt file
    // std::vector <int> edge;
    meshClass::import_matrix_from_txt_file("edge1Labels.txt",edge1);

    // loading point indices of edge from txt file
    // std::vector <int> edge;
    meshClass::import_matrix_from_txt_file("edge2Labels.txt",edge2);
}

int lengthPoints = getMeshLengths(std::vector<int> sizes){
  meshClass::import_matrix_from_txt_file("sizes.txt", sizes);
  return;
}



} // end namespace std
