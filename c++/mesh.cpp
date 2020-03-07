#include "mesh.hpp"
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

        // Structure triangles
        //  - one dimensional vector containing the indices of the Xpoints and Ypoints vectors. 
        //  - each pair of three elements (index%3 = {0,1,2}) forms a triangle.
        const vector<int> triangles_; 

        // Structure triangles
        //  - one dimensional vector containing the indices of the Xpoints and Ypoints vectors. 
        //  - the indices are given in clockwise direction around the edge.
        const vector<int> edge_;

        //vector<vector<T>> 

    public:
    mesh(vector<float> Xpoints, vector<float> Ypoints, vector<int> triangles, vector<int> edge)
    : Xpoints_(Xpoints)
    , Ypoints_(Ypoints)
    , nb_nodes_(Xpoints_.size())
    , nb_elements_(triangles.size()/3)
    , nb_boundary_nodes_(edge_.size())
    , triangles_(triangles)
    , edge_(edge) 
    {
        assert(Xpoints_.size() == Ypoints_.size());
        assert(triangles_.size() % 3 == 0);
    }

    // Checks if point index is on the edge and returns which edge.
    // INPUT: point index 
    // OUTPUT:  - 0 if the point is not on the edge
    //          - 1 for edge 1
    //          - 2 for edge 2
    const int isOnEdge(int index) {
        if (std::find(edge_.begin(), edge_.end(), index) != edge_.end()){
            if (Xpoints_[index] == 0) {
                return 1;
            }
            else {
                return 2;
            }
        }
        else {
            return 0;
        }
    }

    const int getNbElements() {
        return nb_elements_;
    }

    const void getElement(int elementIndex, vector<int> nodes) {
        // return node indices for given element index

    } 

    const void getNodeCoordinates(int nodeIndex, vector<float> coordinates) {
        coordinates[0] = Xpoints_[nodeIndex];
        coordinates[1] = Ypoints_[nodeIndex];
    }

    const int getNbBoundaryNodes() {
        return nb_boundary_nodes_;
    }

    const void getBoundaryNodes(vector<int> nodes, bool boundaryFlag = false) {
        // return nodes indices boundary nodes in the correct order
        // Gamma_1 if boundaryFlag = true, else Gamma_2
    } 

    const int getNbNodes() {
        return nb_nodes_;
    }


}; // end class mesh


// Appends values from a string into a (row)vector and returns the number of elements in the vector (= #columns).
template<class T>
int ReadNumbers( const string & s, vector <T>& v ) {
    istringstream is( s );
    double n;
    //This loop, intuitively, means "keep reading values from is into n, and as long as a value can be read, continue looping." 
    while( is >> n ) {
        // Adds a new element at the end of the vector, after its current last element.
        v.push_back( n );
    }
    return v.size();
}


// Imports al the numbers from a txt file into a one dimensional vector.
template<class T>
void import_matrix_from_txt_file(const char* filename_X, vector <T>& v, int& rows, int& cols){
    
    ifstream file_X;
    string line;
    
    file_X.open(filename_X);
    if (file_X.is_open()) {
        int i=0;
        getline(file_X, line);
        
        
        cols =ReadNumbers( line, v );
        cout << "cols:" << cols << endl;

        // Loop over al the rows of the matrix in the txt file.
        // eof = "end of file."
        while (!file_X.eof()) {
            i = i+1;
            getline(file_X, line);
            ReadNumbers( line, v );
    }
        
        rows=i;
        cout << "rows :" << rows << endl;
        if(rows >32766) cout<< "N must be smaller than MAX_INT";
        
        file_X.close();
    }

    else{
        cout << "file open failed";
    }
    
    // Print out the imported vector.
    cout << "v:" << endl;
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            cout << v[i*cols+j] << "\t" ;
        }
        cout << endl;
    }
}


// Load mesh from text files.
void loadMesh(vector <float> Xpoints, vector <float> Ypoints, vector <int> triangles, vector <int> edge) {
    // loading X-coordinates from txt file
    // vector <float> Xpoints;
    int xnbpoints=0;
    int xpointRows=0;
    import_matrix_from_txt_file("Xpoints.txt",Xpoints,xpointRows,xnbpoints);

    // loading Y-coordinates from txt file
    // vector <float> Ypoints;
    int ynbpoints=0;
    int ypointRows=0;
    import_matrix_from_txt_file("Ypoints.txt",Ypoints,ypointRows,ynbpoints);
    
    // loading point indices of triangles from txt file
    // vector <int> triangles;
    int triangleRows=0;
    int nbtriangles=0;
    import_matrix_from_txt_file("triangleLabels.txt",triangles,triangleRows,nbtriangles);

    // loading point indices of edge from txt file
    // vector <int> edge;
    int edgeRows=0;
    int nbedge=0;
    import_matrix_from_txt_file("edgeLabels.txt",edge,edgeRows,nbedge);
}


} // end namespace std