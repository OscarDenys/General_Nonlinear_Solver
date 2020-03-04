#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

namespace std{



const class mesh{

    private:
        vector<float> Xpoints;
        vector<float> Ypoints;

        // Structure triangles
        //  - one dimensional vector containing the indices of the Xpoints and Ypoints vectors. 
        //  - each pair of three elements (index%3 = {0,1,2}) forms a triangle.
        vector<int> triangles; 

        // Structure triangles
        //  - one dimensional vector containing the indices of the Xpoints and Ypoints vectors. 
        //  - the indices are given in clockwise direction around the edge.
        vector<int> edge;

        //vector<vector<T>> 

    public:
    mesh(vector<float> Xpoints, vector<float> Ypoints, vector<int> triangles, vector<int> edge){
        triangles = triangles;
        Xpoints = Xpoints;
        Ypoints = Ypoints;
        edge = edge;


    }

    const int getNbElements() {
        
    }

    const void getElement(int elementIndex, vector<int> nodes) {
        // return node indices for given element index

    } 

    const void getNodeCoordinates(int nodeIndex, vector<float> coordinates) {


    }

    const int getNbBoundaryNodes() {
        
    }

    const void getBoundaryNodes(vector<int> nodes, boolean boundaryFlag = false) {
        // return nodes indices boundary nodes in the correct order
        // Gamma_1 if boundaryFlag = true, else Gamma_2
    } 

    const int getNbNodes() {

    }


}; //class mesh



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



}// end namespace pear