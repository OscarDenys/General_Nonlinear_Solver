#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
namespace std{



//https://www.codeproject.com/Questions/210023/How-to-read-matrix-from-text-file-in-Cplusplus
//http://www.cplusplus.com/doc/tutorial/files/

class mesh{

    private:
        vector<double> Xpoints;
        vector<double> Ypoints;
        vector<double> triangles; // index in de puntenlijst
        vector<double> edges;

    public:
    mesh(vector<double> Xpoints, vector<double> Ypoints, vector<double> triangles, vector<double> edges){
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
};


//class element{
//
//    public:
//        int label;
//        double x;
//        double y;
//        
//    point(int label, double x, double y){
//        label = label;
//        x = x;
//        y = y;
//    }

    // functie voor edge() -- bool
               //voor neighbours() -- int vector

//};





}// namespace std