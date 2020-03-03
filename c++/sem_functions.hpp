#ifndef sem_functions_hpp
#define sem_functions_hpp
#include "constants.hpp"
#include "mesh.hpp"
#include <Eigen/SparseCore>
#include <cassert>
#include <cmath>
#include <vector>

typedef Eigen::Triplet<double> Trip;

namespace std {

    // respiration kinetics: R - formula (3)

    void evaluateRespiration(int nodeIndex, std::vector<double> prevSol, double Ru, double Rv) {
        // evaluate the respiratory function (formula (3) of assignement) on the node with given index
        // return the values in the doubles Ru and Rv
        // the argument prevSol contains the solution vector of c_i of the previous iteration
        int M = prevSol.size()/2;
        assert(nodeIndex < M);

        Ru = evaluateRu(prevSol[nodeIndex],prevSol[M+nodeIndex]);
        Rv = evaluateRv(prevSol[nodeIndex],prevSol[M+nodeIndex], Ru);
    }

    void evaluateRespiration(int nodeIndex1, int nodeIndex2, std::vector<double> prevSol, double Ru, double Rv) {
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

    // ------------------------

    // line integral evaluation

    void lineIntegral(mesh & myMesh, vector<Trip> & K, Eigen::VectorXd & f) {
        // create matrix K and vector b, where the line integral
        // equals (K*c + b)

        int M = myMesh.getNbNodes();
        int nbBoundaryNodes = myMesh.getNbBoundaryNodes();
        vector<int> boundaryNodes(nbBoundaryNodes);
        myMesh.getBoundaryNodes(boundaryNodes);

        vector<float> nodeCoords(2);
        int currNode, prevNode, nextNode;
        double r_curr, r_prev, r_next, k_prev, k_next;

        for (int i = 1; i < nbBoundaryNodes-1; i++) {

            // node indices
            currNode = boundaryNodes[i];
            prevNode = boundaryNodes[i-1];
            nextNode = boundaryNodes[i+1];
            
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

            K.push_back(Trip(currNode,prevNode, rho_u* k_prev*(r_curr+r_prev)/12 ));
            K.push_back(Trip(currNode,currNode, 
                            rho_u* (k_prev*(3*r_curr+r_prev)/12 + k_next*(3*r_curr+r_next)/12) ) );
            K.push_back(Trip(currNode,nextNode, rho_u* k_next*(r_curr+r_next)/12 ));

            f(currNode) -= rho_u* C_uamb * (k_prev*(2*r_curr+r_prev)/6 + k_next*(2*r_curr+r_next)/6);

            // node indices for v coefficients 
            currNode += M;
            prevNode += M;
            nextNode += M;

            K.push_back(Trip(currNode,prevNode, rho_v* k_prev*(r_curr+r_prev)/12 ));
            K.push_back(Trip(currNode,currNode, 
                            rho_v* (k_prev*(3*r_curr+r_prev)/12 + k_next*(3*r_curr+r_next)/12) ) );
            K.push_back(Trip(currNode,nextNode, rho_v* k_next*(r_curr+r_next)/12 ));

            f(currNode) -= rho_v* C_vamb * (k_prev*(2*r_curr+r_prev)/6 + k_next*(2*r_curr+r_next)/6);

        }
    }

} // namespace std


#endif	