#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "mesh.hpp"
//#ifndef "mesh.hpp"
//#define "mesh.hpp"

namespace std{

    // constructor mesh object
    mesh::mesh(vector<double> Xpoints, vector<double> Ypoints, vector<double> triangles, vector<double> edges):
    //Xpoints_(0), Ypoints_(0), triangles_(0), edges_(0) doesnt fix class redefinition error
    {
        triangles = triangles;
        Xpoints = Xpoints;
        Ypoints = Ypoints;
        edges = edges;

        std::cout << "Xpoints size  = " << Xpoints.size() << std::endl;
        std::cout << "Ypoints size  = " << Ypoints.size() << std::endl;
        //std::cout << "Ypoints   = " << Ypoints << std::endl;

        //std::vector<double> elements(triangles.size()*3);
        //for (int i=0 ; i<triangles.size(); i++){
        //    if (i%3 == 0){
        //        double x = points(2*triangles(i));
        //        double y = points(2*triangles(i)+1);
        //        double edge = 
        //    }
        //}
    }

}// namespace std

//#endif