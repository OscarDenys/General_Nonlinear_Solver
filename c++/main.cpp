//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/
//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "mesh.hpp"
#include "mesh.cpp"

using namespace std;

// Headers 
int ReadNumbers( const string & s, vector <double> & v );
void import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols);

int main(){
    // loading X-coordinates from txt file
    vector <double> Xpoints;
    int xnbpoints=0;
    int xpointRows=0;
    import_matrix_from_txt_file("Xpoints.txt",Xpoints,xpointRows,xnbpoints);

    // loading Y-coordinates from txt file
    vector <double> Ypoints;
    int ynbpoints=0;
    int ypointRows=0;
    import_matrix_from_txt_file("Ypoints.txt",Ypoints,ypointRows,ynbpoints);
    
    // loading point indices of triangles from txt file
    vector <double> triangles;
    int triangleRows=0;
    int nbtriangles=0;
    import_matrix_from_txt_file("triangleLabels.txt",triangles,triangleRows,nbtriangles);

    // loading point indices of edge from txt file
    vector <double> edges;
    int edgeRows=0;
    int nbedge=0;
    import_matrix_from_txt_file("edgeLabels.txt",edges,edgeRows,nbedge);

    std::mesh firstmesh(Xpoints, Ypoints, triangles, edges);



}


// Appends values from a string into a (row)vector and returns the number of elements in the vector (= #columns).
int ReadNumbers( const string & s, vector <double> & v ) {
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
void import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols){
    
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

