#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

namespace std{



const class mesh{

    private:
        vector<double> Xpoints;
        vector<double> Ypoints;

        // Structure triangles & edges:
        //  one dimensional vector containing the indices of the points. 
        //  for triangles each pair of three elements (index%3 = {0,1,2}) forms a triangle.
        vector<int> triangles; 
        vector<int> edge;

    public:
    mesh(vector<double> Xpoints, vector<double> Ypoints, vector<int> triangles, vector<int> edge){
        triangles = triangles;
        Xpoints = Xpoints;
        Ypoints = Ypoints;
        edge = edge;

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

}// namespace std