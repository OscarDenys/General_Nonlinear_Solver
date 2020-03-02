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
        vector<double> triangles; // index in de puntenlijst
        vector<double> edges;

    public:
    mesh(vector<double> Xpoints, vector<double> Ypoints, vector<double> triangles, vector<double> edges);
    
}; //class mesh

}// namespace std