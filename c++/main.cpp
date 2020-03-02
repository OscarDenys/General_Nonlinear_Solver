//http://youngmok.com/c-code-for-reading-unknown-size-matrix-from-text-file/
//http://www.cplusplus.com/doc/tutorial/files/
//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "mesh.hpp"

using namespace std;


int main(){
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

    std::mesh firstmesh(Xpoints, Ypoints, triangles, edges);


}
