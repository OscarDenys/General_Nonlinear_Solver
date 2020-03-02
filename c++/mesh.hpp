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

        //std::cout << "Xpoints size  = " << Xpoints.size() << std::endl;
        //std::cout << "Ypoints size  = " << Ypoints.size() << std::endl;
        

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



}// namespace std