#include <cassert>
#include <cmath>
#include <vector>
#include "constants.hpp"
#include "sibren_functions.hpp"

namespace std {

void createLinearSystem(mesh originalMesh, std::vector<double> K, std::vector<double> f){

  // Assert all elements in given K and f are zero.


  // nodes --> list met nodes van mesh (label, x, y)
  // elements --> list met elements van mesh (node1, node2, node3, boundary(boolean))


  std::vector<int> currentElement(3);
  std::vector<float> P1(2);
  std::vector<float> P2(2);
  std::vector<float> P3(2);
  int n1;
  int n2;
  int n3;
  double det_jac;
  double temp;
  std::vector<double> commonPart(2);
  double resultU;
  double resultV;
  int nbNodes = originalMesh.getNbNodes(); // TODO check if this works;

  for (int element = 0; element < originalMesh.getNbElements(); element++){
    originalMesh.getElement(element, currentElement);
    n1 = currentElement[0];
    n2 = currentElement[1];
    n3 = currentElement[2];
    originalMesh.getNodeCoordinates(n1, P1);
    originalMesh.getNodeCoordinates(n2, P2);
    originalMesh.getNodeCoordinates(n3, P3);

    det_jac = P1[1]*P2[0] - P1[0]*P2[1] + P1[0]*P3[1] - P1[1]*P3[0] - P2[0]*P3[1] + P2[1]*P3[0]; // Bewerking op single prec floats opslaan in double precision float mag???
    temp = (P1[0] + P2[0] + P3[0])/(6*det_jac);

    // ---------------------------------------------------------------------------------
    // This code snippet processes both the first integral of (5) and (6)
    // TODO: Integraal 1; eventueel deelfunctie maken hiervoor???
    // TODO: herschrijf deze assignments zodra we weten hoe sparse matrix voorgesteld wordt.
    // TODO: berekening commonPart vector kan geschreven worden in 1 lijn met vector bewerkingen?
    //    door symmetrie van bewerkingen in r en z....


    // Node 1
    commonPart[0] = pow(P3[0]-P2[0],2);
    commonPart[1] = pow(P3[1]-P2[1],2);
    // commonPart is a vector which contains the part to multiply with sigma_r in [0] and
    // the part with sigma_z in [1]...
    applySigmaAndAddCommonPart(result, commonPart);
    // result is vector with [result_u, result_v]
    K[n1-1,n1-1] = K[n1-1,n1-1] + temp*result[0];
    K[nbNodes+n1-1,nbNodes+n1-1] = K[nbNodes+n1-1,nbNodes+n1-1] + temp*result[1];


    commonPart[0] = -(P3[0]-P2[0])*(P3[0]-P1[0]);
    commonPart[1] = (P2[1]-P3[1])*(P3[1]-P1[1]);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n1-1,n2-1] = K[n1-1,n2-1] + temp*result[0];
    K[nbNodes+n1-1,nbNodes+n2-1] = K[nbNodes+n1-1,nbNodes+n2-1] + temp*result[1];

    commonPart[0] = (P3[0]-P2[0])*(P2[0]-P1[0]);
    commonPart[1] = -(P2[1]-P3[1])*(P2[1]-P1[1]);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n1-1,n3-1] = K[n1-1,n3-1] + temp*result[0];
    K[nbNodes+n1-1,nbNodes+n3-1] = K[nbNodes+n1-1,nbNodes+n3-1] + temp*result[1];



    // Node 2
    commonPart[0] = -(P3[0]-P2[0])*(P3[0]-P1[0]);
    commonPart[1] = (P2[1]-P3[1])*(P3[1]-P1[1]);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n2-1,n1-1] = K[n2-1,n1-1] + temp*result[0];
    K[nbNodes+n2-1,nbNodes+n1-1] = K[nbNodes+n2-1,nbNodes+n1-1] + temp*result[1];

    commonPart[0] = pow(P3[0]-P1[0], 2);
    commonPart[1] = pow(P3[1]-P1[1], 2);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n2-1,n2-1] = K[n2-1,n2-1] + temp*result[0];
    K[nbNodes+n2-1,nbNodes+n2-1] = K[nbNodes+n2-1,nbNodes+n2-1] + temp*result[1];

    commonPart[0] = -(P3[0]-P1[0])*(P2[0]-P1[0]);
    commonPart[1] = -(P3[1]-P1[1])*(P2[1]-P1[1]);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n2-1,n3-1] = K[n2-1,n3-1] + temp*result[0];
    K[nbNodes+n2-1,nbNodes+n3-1] = K[nbNodes+n2-1,nbNodes+n3-1] + temp*result[1];


      // Node 3
    commonPart[0] = (P3[0]-P2[0])*(P2[0]-P1[0]);
    commonPart[1] = -(P2[1]-P3[1])*(P2[1]-P1[1]);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n3-1,n1-1] = K[n3-1,n1-1] + temp*result[0];
    K[nbNodes+n3-1,nbNodes+n1-1] = K[nbNodes+n3-1,nbNodes+n1-1] + temp*result[1];

    commonPart[0] = -(P3[0]-P1[0])*(P2[0]-P1[0]);
    commonPart[1] = -(P3[1]-P1[1])*(P2[1]-P1[1]);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n3-1,n2-1] = K[n3-1,n2-1] + temp*result[0];
    K[nbNodes+n3-1,nbNodes+n2-1] = K[nbNodes+n3-1,nbNodes+n2-1] + temp*result[1];

    commonPart[0] = pow(P2[0]-P1[0], 2);
    commonPart[1] = pow(P2[1]-P1[1], 2);
    applySigmaAndAddCommonPart(result, commonPart);
    K[n3-1,n3-1] = K[n3-1,n3-1] + temp*result[0];
    K[nbNodes+n3-1,nbNodes+n3-1] = K[nbNodes+n3-1,nbNodes+n3-1] + temp*result[1];


  }


}

void applySigmaAndAddCommonPart(std::vector<double> result, std::vector<double> commonPart){
  // TODO: volgens wiskundige afleiding moet sigma_uz maal r gedeelte en andersom, dit is niet hoe wij het in MATLAB doen!!!
  result[0] = sigma_uz*commonPart[0] + sigma_ur*commonPart[1];
  result[1] = sigma_vz*commonPart[0] + sigma_vr*commonPart[1];

  return;
}



}
