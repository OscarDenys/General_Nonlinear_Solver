#include <iostream>
#include <fstream>
#include <string>
using namespace std;
    
    
    
    
    
class mesh {
    private:
        vector<double> points;
        vector<double> triangles; // labels zijn in volgorde van de punten lijst, dus 1 is het eerste punt
        vector<double> edges;

    public:
    mesh(vector<double> points, vector<double> triangles, vector<double> edges);

}; // class mesh