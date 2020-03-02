#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

namespace std{



class mesh{

    private:
        vector<double> Xpoints;
        vector<double> Ypoints;

        // Structure triangles & edges:
        //  one dimensional vector containing the indices of the points. 
        //  for triangles each pair of three elements (index%3 = {0,1,2}) forms a triangle.
        vector<double> triangles; 
        vector<double> edge;

    public:
    mesh(vector<double> Xpoints, vector<double> Ypoints, vector<double> triangles, vector<double> edge){
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
    
}; //class mesh

}// namespace std