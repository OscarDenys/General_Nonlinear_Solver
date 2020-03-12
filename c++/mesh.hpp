#ifndef mesh_hpp
#define mesh_hpp
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cassert>

namespace meshClass{

class mesh{

    private:
      const std::vector<float> Xpoint_;
      const std::vector<float> Ypoint_;
      std::vector<int> triangles_;
      std::vector<int> edge1_;
      std::vector<int> edge2_;

    public:
    mesh(std::vector<float> Xpoint, std::vector<float> Ypoint, std::vector<int> triangles, std::vector<int> edge1, std::vector<int> edge2);


    const int getNbElements();
    const int getNbNodes();
    const void getElement(int elementIndex, std::vector<int> &nodeIndices);
    const void getNodeCoordinates(int nodeIndex, std::vector<float> coordinates);
    const int getNbBoundaryNodes(bool edge2 = true);
    const void getBoundaryNodes(std::vector<int>& nodeIndices, bool edge2 = true);


}; // end class mesh


// Appends values from a string into a (row)vector and returns the number of elements in the vector (= #columns).
template<class T>
int ReadNumbers( const std::string & s, std::vector <T>& v );

// Imports al the numbers from a txt file into a one dimensional vector.
template<class T>
void import_matrix_from_txt_file(const char* filename_X, std::vector<T>& v);

// Load mesh from text files.
void loadMesh(std::vector<float> Xpoint, std::vector<float> Ypoint, std::vector<int> triangles, std::vector<int> edge1, std::vector<int> edge2);

int lengthPoints = getMeshLengths(std::vector<int> sizes);

}


#endif
