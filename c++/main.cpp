//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/


#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "mesh.hpp"
#include "mesh.cpp"

using namespace std;

int ReadNumbers( const string & s, vector <double> & v );
void import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols);

int main(){
    vector <double> Xpoints;
    int xnbpoints=0;
    int xpointRows=0;
    import_matrix_from_txt_file("Xpoints.txt",Xpoints,xpointRows,xnbpoints);

    vector <double> Ypoints;
    int ynbpoints=0;
    int ypointRows=0;
    import_matrix_from_txt_file("Ypoints.txt",Ypoints,ypointRows,ynbpoints);
    //std::cout << "points size  = " << points.size() << std::endl;

    vector <double> triangles;
    int triangleRows=0;
    int nbtriangles=0;
    import_matrix_from_txt_file("triangleLabels.txt",triangles,triangleRows,nbtriangles);

    vector <double> edges;
    int edgeRows=0;
    int nbedge=0;
    import_matrix_from_txt_file("edgeLabels.txt",edges,edgeRows,nbedge);

    std::mesh firstmesh(Xpoints, Ypoints, triangles, edges);



}

int ReadNumbers( const string & s, vector <double> & v ) {
    istringstream is( s );
    double n;
    while( is >> n ) {
        v.push_back( n );
    }
    return v.size();
}



void import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols){
    
    ifstream file_X;
    string line;
    
    file_X.open(filename_X);
    if (file_X.is_open()) {
        int i=0;
        getline(file_X, line);
        
        
        cols =ReadNumbers( line, v );
        cout << "cols:" << cols << endl;
        
       // while ( getline (file_X,line) ){
        //    ReadNumbers( line, v );
       // }

        while (!file_X.eof()) {
            i = i+1;
            getline(file_X, line);
            ReadNumbers( line, v );
    }

        //for ( i=1;i<32767;i++){
         //   getline(file_X, line);
         //   if ( line == 0 ) break;
         //   ReadNumbers( line, v );
       // }
        
        rows=i;
        cout << "rows :" << rows << endl;
        if(rows >32766) cout<< "N must be smaller than MAX_INT";
        
        file_X.close();
    }

    else{
        cout << "file open failed";
    }
    
    cout << "v:" << endl;
    //cout << v.size(1) << endl;//
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            cout << v[i*cols+j] << "\t" ;
        }
        cout << endl;
    }
}

//}//namespace pme