#ifndef mesh_hpp
#define mesh_hpp
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cassert>

namespace std{



class mesh{

    private:
        const vector<float> Xpoints_;
        const vector<float> Ypoints_;
        const int nb_nodes_;
        const int nb_elements_;
        const int nb_boundary_nodes_;
        const vector<int> triangles_; 
        const vector<int> edge_;

    public:
    mesh(vector<float> Xpoints, vector<float> Ypoints, vector<int> triangles, vector<int> edge);

    // Checks if point index is on the edge and returns which edge.
    // INPUT: point index 
    // OUTPUT:  - 0 if the point is not on the edge
    //          - 1 for edge 1
    //          - 2 for edge 2
    const int isOnEdge(int index);

    const int getNbElements();

    const void getElement(int elementIndex, vector<int> nodes);

    const void getNodeCoordinates(int nodeIndex, vector<float> coordinates);

    const int getNbBoundaryNodes();

    const void getBoundaryNodes(vector<int> nodes, bool boundaryFlag = false);

    const int getNbNodes();

}; // end class mesh


// Appends values from a string into a (row)vector and returns the number of elements in the vector (= #columns).
template<class T>
int ReadNumbers( const string & s, vector <T>& v );

// Imports al the numbers from a txt file into a one dimensional vector.
template<class T>
void import_matrix_from_txt_file(const char* filename_X, vector <T>& v, int& rows, int& cols);

// Load mesh from text files.
void loadMesh(vector <float> Xpoints, vector <float> Ypoints, vector <int> triangles, vector <int> edge);


} 


#endif