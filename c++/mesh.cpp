//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/
//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus


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
        const vector<float> Xpoint_;
        const vector<float> Ypoint_;
        const int nb_nodes_;
        const int nb_elements_;
        const int nb_boundary_nodes_;

        // Structure triangles
        //  - one dimensional vector containing the indices of the Xpoints and Ypoints vectors. 
        //  - each pair of three elements (index%3 = {0,1,2}) forms a triangle.
        const vector<int> triangles_; 

        // Structure edge
        //  - one dimensional vector containing the indices of the Xpoint and Ypoint vectors. 
        //  - the indices are given in clockwise direction around the edge (starting from the top of the pear).
        const vector<int> edge_;
        const vector<int> edge1_;
        const vector<int> edge2_;

        const vector<vector<int>> element_;

    public:
    mesh(vector<float> Xpoint, vector<float> Ypoint, vector<int> triangles, vector<int> edge)
    : Xpoint_(Xpoint)
    , Ypoint_(Ypoint)
    , nb_nodes_(Xpoint_.size())
    , nb_elements_(triangles.size()/3)
    , nb_boundary_nodes_(edge_.size())
    , triangles_(triangles)
    , edge_(edge) 
    {
        assert(Xpoint_.size() == Ypoint_.size());
        assert(triangles_.size() % 3 == 0);

        const vector<vector<int>> element_(nb_elements_);
        for (int i = 0; i < triangles.size(); i++){
            if (i%3 ==0){
                element_(i/3).push_back(triangles_(i));
                element_(i/3).push_back(triangles_(i+1));
                element_(i/3).push_back(triangles_(i+2));
            }
        }

        const vector<int> edge1_;
        const vector<int> edge2_;
        for (int i = 0; i < edge_.size(); i++){
            if ( Xpoint_(edge_(i)) == 0 ) {
                edge1_.push_back(edge(i));
            }
            else {
                edge2_.push_back(edge(i));
            }
        }

    }

    // Checks if point index is on the edge and returns which edge.
    // INPUT: point index 
    // OUTPUT:  - 0 if the point is not on the edge
    //          - 1 for edge 1
    //          - 2 for edge 2
    const int isOnEdge(int index) {
        if (std::find(edge_.begin(), edge_.end(), index) != edge_.end()){
            if (Xpoint_[index] == 0) {
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

    const void getElement(int elementIndex, vector<int> nodeIndices) {
        // return node indices for given element index
        assert (nodeIndices.size() == 3);
        nodeIndices.swap(element_(elementIndex));
    } 

    const void getNodeCoordinates(int nodeIndex, vector<float> coordinates) {
        coordinates[0] = Xpoint_[nodeIndex];
        coordinates[1] = Ypoint_[nodeIndex];
    }

    // Returns total number of boundary nodes without argument. 
    // Returns number of nodes on edge 1 or 2 if specified.
    const int getNbBoundaryNodes(int edge1or2 = 0) {
        if (edge1or2 ==0){
            return nb_boundary_nodes_;
        }
        else if (edge1or2 ==1){
            return edge1_.size();
        }
        else if (edge1or2 ==2){
            return edge2_.size();
        }
    }

    const void getBoundaryNodes(vector<int> nodeIndices, bool firstBoundary = false) {
        // return nodes indices boundary nodes in the correct order
        // Gamma_1 if firstBoundary = true, else Gamma_2
        if (firstBoundary){
            nodeIndices.swap(edge1_);
        }
        else {
            nodeIndices.swap(edge2_);
        }
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
void loadMesh(vector <float> Xpoint, vector <float> Ypoint, vector <int> triangles, vector <int> edge) {
    // loading X-coordinates from txt file
    // vector <float> Xpoints;
    int xnbpoints=0;
    int xpointRows=0;
    import_matrix_from_txt_file("Xpoints.txt",Xpoint,xpointRows,xnbpoints);

    // loading Y-coordinates from txt file
    // vector <float> Ypoints;
    int ynbpoints=0;
    int ypointRows=0;
    import_matrix_from_txt_file("Ypoints.txt",Ypoint,ypointRows,ynbpoints);
    
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