#include "sem_functions.hpp"
#include "constants.hpp"
#include "mesh.hpp"
#include <Eigen/SparseCore>
#include <cassert>
#include <cmath>
#include <vector>



typedef Eigen::Triplet<double> Trip;

namespace std {

    // determinant of the Jacobian
    //
    double detJac(vector<float> P1, vector<float> P2, vector<float> P3) {
        return double(P1[1]*P2[0] - P1[0]*P2[1] + P1[0]*P3[1] - P1[1]*P3[0] - P2[0]*P3[1] + P2[1]*P3[0]); // Bewerking op single prec floats opslaan in double precision float mag???
    }


    // respiration kinetics: R - formula (3)
    //
    void evaluateRespiration(int nodeIndex, Eigen::VectorXd &prevSol, double &Ru, double &Rv) {
        // evaluate the respiratory function (formula (3) of assignement) on the node with given index
        // return the values in the doubles Ru and Rv
        // the argument prevSol contains the solution vector of c_i of the previous iteration
        int M = prevSol.size()/2;
        assert(nodeIndex < M);

        Ru = evaluateRu(prevSol[nodeIndex],prevSol[M+nodeIndex]);
        Rv = evaluateRv(prevSol[nodeIndex],prevSol[M+nodeIndex], Ru);
    }

    void evaluateRespiration(int nodeIndex1, int nodeIndex2, Eigen::VectorXd &prevSol, double &Ru, double &Rv) {
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
    //
    // ---------------------------------------

    // second (linearised) integral evaluation
    //
    void integral2lin(mesh & myMesh, vector<Trip> & K_lin, Eigen::VectorXd &f_lin) {
        // create matrix K_lin and vector f_lin, where the second integral approximates
        // (K_lin*c + f_lin)
        //
        //      K_lin not used in this linearisation
        //

        double Ru = evaluateRu(C_uamb, C_vamb);
        double Rv = evaluateRv(C_uamb, C_vamb, Ru);
        int M = myMesh.getNbNodes();
        double det_jac;
        int n1, n2, n3;
        vector<int> currElem(3);
        vector<float> P1(2), P2(2), P3(3);

        int nbElements = myMesh.getNbElements();

        for (int elem = 0; elem < nbElements; elem++) {
            myMesh.getElement(elem, currElem);
            n1 = currElem[0];
            n2 = currElem[1];
            n3 = currElem[2];
            myMesh.getNodeCoordinates(n1, P1);
            myMesh.getNodeCoordinates(n2, P2);
            myMesh.getNodeCoordinates(n3, P3);

            det_jac = detJac(P1,P2,P3);

            f_lin[n1]   -= det_jac * Ru * (2*P1[0] + P2[0] + P3[0]) / 24;
            f_lin[n2]   -= det_jac * Ru * (P1[0] + 2*P2[0] + P3[0]) / 24;
            f_lin[n3]   -= det_jac * Ru * (P1[0] + P2[0] + 2*P3[0]) / 24;
            f_lin[n1+M] -= det_jac * Rv * (2*P1[0] + P2[0] + P3[0]) / 24;
            f_lin[n2+M] -= det_jac * Rv * (P1[0] + 2*P2[0] + P3[0]) / 24;
            f_lin[n3+M] -= det_jac * Rv * (P1[0] + P2[0] + 2*P3[0]) / 24;
        }

    }
    //
    // ------------------------

    // line integral evaluation
    //
    void integral3(mesh & myMesh, vector<Trip> & K, Eigen::VectorXd & f) {
        // create matrix K and vector f, where the line integral equals
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

    }
    //
    // ------------------------

} // namespace std
