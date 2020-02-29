#include <cassert>
#include <cmath>
#include <vector>
#include "sibren_functions.hpp"

namespace pear{

void createLinearSystem(mesh originalMesh, std::vector<double> K, std::vector<double> f, std::vector<float> setup){
  // Assert all elements in given K and f are zero.

  // nodes --> list met nodes van mesh (label, x, y)
  // elements --> list met elements van mesh (node1, node2, node3, boundary(boolean))


  std::vector<int> currentElement(3);
  boolean isBoundary;
  std::vector<float> P1(2);
  std::vector<float> P2(2);
  std::vector<float> P3(2);
  double det_jac;
  double temp;
  double result;
  for (int element = 0; element < elements.size(); element++){
    currentElement = elements[element];
    isBoundary = currentElement[3];
    P1 = nodes[currentElement[1]];
    P2 = nodes[currentElement[2]];
    P3 = nodes[currentElement[3]];

    det_jac = P1[2]*P2[1] - P1[1]*P2[2] + P1[1]*P3[2] - P1[2]*P3[1] - P2[1]*P3[2] + P2[2]*P3[1]; // Bewerking op single prec floats opslaan in double precision float mag???
    temp = (P1[1] + P2[1] + P3[1])/(6*det_jac);

    // TODO: Integraal 1; eventueel deelfunctie maken hiervoor???
    // TODO: herschrijf deze assignments zodra we weten hoe sparse matrix voorgesteld wordt.
      // Node 1
    result = temp * (sigma_ur*pow(P3[1]-P2[1],2) + sigma_uz*pow(P3[2]-P2[2],2)) ;
    K[currentElement[1],currentElement[1]] = K[currentElement[1],currentElement[1]] + result;
    result = temp * (-sigma_ur*(P3[1]-P2[1])(P3[1]-P1[1]) + sigma_uz*(P2[2]-P3[2])*(P3[2]-P1[2])) ;
    K[currentElement[1],currentElement[2]] = K[currentElement[1],currentElement[2]] + result;
    result = temp * (sigma_ur*(P3[1]-P2[1])(P2[1]-P1[1]) - sigma_uz*(P2[2]-P3[2])*(P2[2]-P1[2])) ;
    K[currentElement[1],currentElement[3]] = K[currentElement[1],currentElement[3]] + result;

      // Node 2
    result = temp * (-sigma_ur*(P3[1]-P2[1])(P3[1]-P1[1]) + sigma_uz*(P2[2]-P3[2])*(P3[2]-P1[2])) ;
    K[currentElement[2],currentElement[1]] = K[currentElement[2],currentElement[1]] + result;
    result = temp * (sigma_ur*pow(P3[1]-P1[1],2) + sigma_uz*pow(P3[2]-P1[2],2)) ;
    K[currentElement[2],currentElement[2]] = K[currentElement[2],currentElement[2]] + result;
    result = temp * (-sigma_ur*(P3[1]-P1[1])(P2[1]-P1[1]) - sigma_uz*(P3[2]-P1[2])*(P2[2]-P1[2])) ;
    K[currentElement[2],currentElement[3]] = K[currentElement[2],currentElement[3]] + result;

      // Node 3
    result = temp * (sigma_ur*(P3[1]-P2[1])(P2[1]-P1[1]) - sigma_uz*(P2[2]-P3[2])*(P2[2]-P1[2])) ;
    K[currentElement[2],currentElement[1]] = K[currentElement[2],currentElement[1]] + result;
    result = temp *  (-sigma_ur*(P3[1]-P1[1])(P2[1]-P1[1]) - sigma_uz*(P3[2]-P1[2])*(P2[2]-P1[2])) ;
    K[currentElement[2],currentElement[2]] = K[currentElement[2],currentElement[2]] + result;
    result = temp * (sigma_ur*pow(P2[1]-P1[1],2) + sigma_uz*pow(P2[2]-P1[2],2)) ;
    K[currentElement[2],currentElement[3]] = K[currentElement[2],currentElement[3]] + result;



  }
}



}
