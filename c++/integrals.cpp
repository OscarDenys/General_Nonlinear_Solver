#include <cassert>
#include <cmath>
#include <vector>
#include "integrals.hpp"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "mesh.hpp"
#include "constants.hpp"


typedef Eigen::Triplet<double> Trip;

// typedef to switch efficiently between float and double precision would be a good idea.
//  Because from the moment a number of single precision is used in a calculation this imediately makes the result single precision.
typedef double prec;
//typedef float prec;

namespace std {

    void integral1(std::mesh &myMesh, std::vector<Trip> &K){
        
        // Preallocate arrays in memory: 
        std::vector<int> currentElement(3);
        std::vector<float> P1(2), P2(2), P3(2);
        int n1, n2, n3;
        double det_jac, temp;
        std::vector<double> commonPart(2), result(2);
        int nbNodes = myMesh.getNbNodes(); 

        // Iterate over all elements and build stiffness matrix: 
        for (int element = 0; element < myMesh.getNbElements(); element++){
        myMesh.getElement(element, currentElement);
        n1 = currentElement[0];
        n2 = currentElement[1];
        n3 = currentElement[2];
        myMesh.getNodeCoordinates(n1, P1);
        myMesh.getNodeCoordinates(n2, P2);
        myMesh.getNodeCoordinates(n3, P3);

        det_jac = detJac(P1,P2,P3);
        
        temp = (P1[0] + P2[0] + P3[0])/(6*det_jac);

        
        // Node 1
        // commonPart is a vector which contains the part to multiply with sigma_r in [0] and
        // the part with sigma_z in [1]...
        commonPart[0] = pow(P3[0]-P2[0],2);
        commonPart[1] = pow(P3[1]-P2[1],2);
        applySigmaAndAddCommonPart(result, commonPart);
        // result is vector with [result_u, result_v]
        K.push_back(Trip(n1, n1, temp*result[0]));
        K.push_back(Trip(nbNodes+n1, nbNodes+n1, temp*result[1]));

        commonPart[0] = -(P3[0]-P2[0])*(P3[0]-P1[0]);
        commonPart[1] = (P2[1]-P3[1])*(P3[1]-P1[1]);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n1, n2, temp*result[0]));
        K.push_back(Trip(nbNodes+n1, nbNodes+n2, temp*result[1]));

        commonPart[0] = (P3[0]-P2[0])*(P2[0]-P1[0]);
        commonPart[1] = -(P2[1]-P3[1])*(P2[1]-P1[1]);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n1, n3, temp*result[0]));
        K.push_back(Trip(nbNodes+n1, nbNodes+n3, temp*result[1]));

        // Node 2
        commonPart[0] = -(P3[0]-P2[0])*(P3[0]-P1[0]);
        commonPart[1] = (P2[1]-P3[1])*(P3[1]-P1[1]);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n2, n1, temp*result[0]));
        K.push_back(Trip(nbNodes+n2, nbNodes+n1, temp*result[1]));

        commonPart[0] = pow(P3[0]-P1[0], 2);
        commonPart[1] = pow(P3[1]-P1[1], 2);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n2, n2, temp*result[0]));
        K.push_back(Trip(nbNodes+n2, nbNodes+n2, temp*result[1]));

        commonPart[0] = -(P3[0]-P1[0])*(P2[0]-P1[0]);
        commonPart[1] = -(P3[1]-P1[1])*(P2[1]-P1[1]);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n2, n3, temp*result[0]));
        K.push_back(Trip(nbNodes+n2, nbNodes+n3, temp*result[1]));

        // Node 3
        commonPart[0] = (P3[0]-P2[0])*(P2[0]-P1[0]);
        commonPart[1] = -(P2[1]-P3[1])*(P2[1]-P1[1]);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n3, n1, temp*result[0]));
        K.push_back(Trip(nbNodes+n3, nbNodes+n1, temp*result[1]));

        commonPart[0] = -(P3[0]-P1[0])*(P2[0]-P1[0]);
        commonPart[1] = -(P3[1]-P1[1])*(P2[1]-P1[1]);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n3, n2, temp*result[0]));
        K.push_back(Trip(nbNodes+n3, nbNodes+n2, temp*result[1]));

        commonPart[0] = pow(P2[0]-P1[0], 2);
        commonPart[1] = pow(P2[1]-P1[1], 2);
        applySigmaAndAddCommonPart(result, commonPart);
        K.push_back(Trip(n3, n3, temp*result[0]));
        K.push_back(Trip(nbNodes+n3, nbNodes+n3, temp*result[1]));
        } // integral1()
    }

    void integral2lin(mesh &myMesh, vector<Trip> & K_lin, Eigen::VectorXd &f_lin) {
        // Linearization of integral 2: H(C_amb)

        // Preallocate arrays in memory:
        double Ru = evaluateRu(C_uamb, C_vamb);
        double Rv = evaluateRv(C_uamb, C_vamb, Ru);
        int M = myMesh.getNbNodes();
        double det_jac;
        int n1, n2, n3;
        vector<int> currElem(3);
        vector<float> P1(2), P2(2), P3(2);

        int nbElements = myMesh.getNbElements();

        // Iterate over elements and fill load vector:
        for (int elem = 0; elem < nbElements; elem++) {
            myMesh.getElement(elem, currElem);
            n1 = currElem[0];
            n2 = currElem[1];
            n3 = currElem[2];
            myMesh.getNodeCoordinates(n1, P1);
            myMesh.getNodeCoordinates(n2, P2);
            myMesh.getNodeCoordinates(n3, P3);

            det_jac = detJac(P1,P2,P3);

            f_lin[n1]   += det_jac * Ru * (2*P1[0] + P2[0] + P3[0]) / 24;
            f_lin[n2]   += det_jac * Ru * (P1[0] + 2*P2[0] + P3[0]) / 24;
            f_lin[n3]   += det_jac * Ru * (P1[0] + P2[0] + 2*P3[0]) / 24;
            f_lin[n1+M] -= det_jac * Rv * (2*P1[0] + P2[0] + P3[0]) / 24;
            f_lin[n2+M] -= det_jac * Rv * (P1[0] + 2*P2[0] + P3[0]) / 24;
            f_lin[n3+M] -= det_jac * Rv * (P1[0] + P2[0] + 2*P3[0]) / 24;
        };
    } // integral2lin()

    Eigen::ArrayXd integral2nonlinear(Eigen::ArrayXd &C, std::mesh &myMesh) {
        // Evaluation of nonlinear H(C):

        // Preallocate arrays in memory: 
        Eigen::ArrayXd H(C.size());  
        double Ru12, Ru13, Ru23;
        double Rv12, Rv13, Rv23;
        int M = myMesh.getNbNodes();
        double det_jac;
        int n1, n2, n3;
        vector<int> currElem(3);
        vector<float> P1(2), P2(2), P3(2);

        int nbElements = myMesh.getNbElements();

        // Iterate over all elements and fill load vector: 
        for (int elem = 0; elem < nbElements; elem++) {
            myMesh.getElement(elem, currElem);
            n1 = currElem[0];
            n2 = currElem[1];
            n3 = currElem[2];
            
            myMesh.getNodeCoordinates(n1, P1);
            myMesh.getNodeCoordinates(n2, P2);
            myMesh.getNodeCoordinates(n3, P3);

            det_jac = detJac(P1,P2,P3);

            evaluateRespiration(n1, n2, C, Ru12, Rv12);
            evaluateRespiration(n1, n3, C, Ru13, Rv13);
            evaluateRespiration(n2, n3, C, Ru23, Rv23);
            
            H[n1]   += det_jac * ((P1[0]+P2[0])*Ru12 + (P1[0] + P3[0])*Ru13) / 24;
            H[n2]   += det_jac * ((P1[0]+P2[0])*Ru12 + (P2[0] + P3[0])*Ru23) / 24;
            H[n3]   += det_jac * ((P1[0]+P3[0])*Ru13 + (P2[0] + P3[0])*Ru23) / 24;
            H[n1+M] -= det_jac * ((P1[0]+P2[0])*Rv12 + (P1[0] + P3[0])*Rv13) / 24;
            H[n2+M] -= det_jac * ((P1[0]+P2[0])*Rv12 + (P2[0] + P3[0])*Rv23) / 24;
            H[n3+M] -= det_jac * ((P1[0]+P3[0])*Rv13 + (P2[0] + P3[0])*Rv23) / 24;
        }
        return H;
    } // integral2nonLinear()


    void integral3(mesh &myMesh, std::vector<Trip> & K, Eigen::VectorXd & f) {
        // create matrix K_3 and vector f, where the line integral equals
        // (K*c + f)

        // init
        int M = myMesh.getNbNodes();
        int nbBoundaryNodes = myMesh.getNbBoundaryNodes();
        vector<int> boundaryNodes(nbBoundaryNodes);
        myMesh.getBoundaryNodes(boundaryNodes);
        vector<float> nodeCoords(2);
        int currNode, prevNode, nextNode;
        double r_curr, r_prev, r_next, k_prev, k_next;

        // only evaluate section (1.6.2) from mathematical derivation for first node
        currNode = boundaryNodes[0];
        nextNode = boundaryNodes[1];

        // value initialisation
        myMesh.getNodeCoordinates(nextNode, nodeCoords);
        r_next = nodeCoords[0];
        k_next = nodeCoords[1];     // temporarily store z-value here
        myMesh.getNodeCoordinates(currNode, nodeCoords);
        r_curr = nodeCoords[0];
        k_next = sqrt( pow(r_next-r_curr,2) + pow(k_next-nodeCoords[1],2) );

        // storing the u-coefficients in the matrix
        K.push_back(Trip(currNode,currNode,
                        rho_u* k_next*(3*r_curr+r_next)/12) );
        K.push_back(Trip(currNode,nextNode, rho_u* k_next*(r_curr+r_next)/12 ));
        f(currNode) -= rho_u* C_uamb * k_next*(2*r_curr+r_next)/6;

        // node indices for v coefficients
        currNode += M;
        nextNode += M;
        // storing the v-coefficients in the matrix
        K.push_back(Trip(currNode,currNode,
                        rho_v* k_next*(3*r_curr+r_next)/12) );
        K.push_back(Trip(currNode,nextNode, rho_u* k_next*(r_curr+r_next)/12 ));
        f(currNode) -= rho_v* C_vamb * k_next*(2*r_curr+r_next)/6;


        // evaluate section (1.6.30) for all but the first and last boundary node
        for (int i = 1; i < nbBoundaryNodes-1; i++) {
            // TODO: think about how this code could be optimised, because every node will be fetched
            //      three times. Once as prevNode, once as currNode, once as nextNode.

            // node indices
            currNode = boundaryNodes[i];
            prevNode = boundaryNodes[i-1];
            nextNode = boundaryNodes[i+1];

            // value initialisation
            myMesh.getNodeCoordinates(prevNode, nodeCoords);
            r_prev = nodeCoords[0];
            k_prev = nodeCoords[1];     // temporarily store z-value here
            myMesh.getNodeCoordinates(nextNode, nodeCoords);
            r_next = nodeCoords[0];
            k_next = nodeCoords[1];     // temporarily store z-value here
            myMesh.getNodeCoordinates(currNode, nodeCoords);
            r_curr = nodeCoords[0];
            k_prev = sqrt( pow(r_prev-r_curr,2) + pow(k_prev-nodeCoords[1],2) );
            k_next = sqrt( pow(r_next-r_curr,2) + pow(k_next-nodeCoords[1],2) );

            // storing the u-coefficients in the matrix
            K.push_back(Trip(currNode,prevNode, rho_u* k_prev*(r_curr+r_prev)/12 ));
            K.push_back(Trip(currNode,currNode,
                            rho_u* (k_prev*(3*r_curr+r_prev)/12 + k_next*(3*r_curr+r_next)/12) ) );
            K.push_back(Trip(currNode,nextNode, rho_u* k_next*(r_curr+r_next)/12 ));
            f(currNode) -= rho_u* C_uamb * (k_prev*(2*r_curr+r_prev)/6 + k_next*(2*r_curr+r_next)/6);

            // node indices for v coefficients
            currNode += M;
            prevNode += M;
            nextNode += M;

            // storing the v-coefficients in the matrix
            K.push_back(Trip(currNode,prevNode, rho_v* k_prev*(r_curr+r_prev)/12 ));
            K.push_back(Trip(currNode,currNode,
                            rho_v* (k_prev*(3*r_curr+r_prev)/12 + k_next*(3*r_curr+r_next)/12) ) );
            K.push_back(Trip(currNode,nextNode, rho_v* k_next*(r_curr+r_next)/12 ));
            f(currNode) -= rho_v* C_vamb * (k_prev*(2*r_curr+r_prev)/6 + k_next*(2*r_curr+r_next)/6);

        }


        // only evaluate section (1.6.1) from mathematical derivation for last node
        currNode = boundaryNodes[nbBoundaryNodes-1];
        prevNode = boundaryNodes[nbBoundaryNodes-2];

        // value initialisation
        myMesh.getNodeCoordinates(prevNode, nodeCoords);
        r_prev = nodeCoords[0];
        k_prev = nodeCoords[1];     // temporarily store z-value here
        myMesh.getNodeCoordinates(currNode, nodeCoords);
        r_curr = nodeCoords[0];
        k_prev = sqrt( pow(r_prev-r_curr,2) + pow(k_prev-nodeCoords[1],2) );

        // storing the u-coefficients in the matrix
        K.push_back(Trip(currNode,currNode,
                        rho_u* k_prev*(3*r_curr+r_prev)/12) );
        K.push_back(Trip(currNode,prevNode, rho_u* k_prev*(r_curr+r_prev)/12 ));
        f(currNode) -= rho_u* C_uamb * k_prev*(2*r_curr+r_prev)/6;

        // node indices for v coefficients
        currNode += M;
        prevNode += M;
        // storing the v-coefficients in the matrix
        K.push_back(Trip(currNode,currNode,
                        rho_v* k_prev*(3*r_curr+r_prev)/12) );
        K.push_back(Trip(currNode,prevNode, rho_u* k_prev*(r_curr+r_prev)/12 ));
        f(currNode) -= rho_v* C_vamb * k_prev*(2*r_curr+r_prev)/6;

    } // integral3()

    double detJac(vector<float> P1, vector<float> P2, vector<float> P3) {
        // Return the determinant of the Jacobian (transformation to local coordinate-system)
        return abs(double(P1[1]*P2[0] - P1[0]*P2[1] + P1[0]*P3[1] - P1[1]*P3[0] - P2[0]*P3[1] + P2[1]*P3[0])); 
    } // detJac                                                                                        


    // Overleaf p3: f1, f2 & f3
    void applySigmaAndAddCommonPart(std::vector<double> &result, std::vector<double> const &commonPart){
        result[0] = sigma_uz*commonPart[0] + sigma_ur*commonPart[1];
        result[1] = sigma_vz*commonPart[0] + sigma_vr*commonPart[1];

        return;
    } // applySigmaAndAddCommonPart()

    // respiration kinetics: R - formula (3)
    void evaluateRespiration(int nodeIndex, Eigen::ArrayXd &prevSol, double &Ru, double &Rv) {
        // evaluate the respiratory function (formula (3) of assignement) on the node with given index
        // return the values in the doubles Ru and Rv
        // the argument prevSol contains the solution vector of c_i of the previous iteration
        int M = prevSol.size()/2;
        assert(nodeIndex < M);

        Ru = evaluateRu(prevSol[nodeIndex],prevSol[M+nodeIndex]);
        Rv = evaluateRv(prevSol[nodeIndex],prevSol[M+nodeIndex], Ru);
    }

    void evaluateRespiration(int nodeIndex1, int nodeIndex2, Eigen::ArrayXd &prevSol, double &Ru, double &Rv) {
        // evaluate the respiratory function (formula (3) of assignement) halfway between two nodes
        // with given index
        // return the values in the doubles Ru and Rv
        // the argument prevSol contains the solution vector of c_i of the previous iteration
        int M = prevSol.size()/2;
        assert(nodeIndex1 < M);
        assert(nodeIndex2 < M);

        double Cu = (prevSol[nodeIndex1]+prevSol[nodeIndex2])/2;
        double Cv = (prevSol[M+nodeIndex1]+prevSol[M+nodeIndex2])/2;

        Ru = evaluateRu(Cu,Cv);
        Rv = evaluateRv(Cu,Cv,Ru);
    }

    double evaluateRu(double Cu, double Cv) {
        // evaluate Ru following formula (3) of the assignement

        return V_mu * Cu / ((K_mu+Cu)*(1+ Cv/K_mv));
    }

    double evaluateRv(double Cu, double Cv, double Ru) {
        // evaluate Rv following formula (3) of the assignement

        return r_q*Ru + V_mfv / (1+ Cu/K_mfu);
    }

    void evaluateCostFunction( Eigen::SparseMatrix<double> &K, Eigen::ArrayXd &f, Eigen::ArrayXd &C, Eigen::ArrayXd &F, mesh &myMesh){
        Eigen::VectorXd Kc = K * C.matrix();
        F = Kc.array() + f + integral2nonlinear(C, myMesh);
    }

} // namespace std
