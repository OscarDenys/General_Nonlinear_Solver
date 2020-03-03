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

        int nbBoundaryNodes = myMesh.getNbBoundaryNodes();
        vector<int> boundaryNodes(nbBoundaryNodes);
        myMesh.getBoundaryNodes(boundaryNodes);

        int currentNode, prevNode, nextNode;

        for (int i = 1; i < nbBoundaryNodes-1; i++) {
            currentNode = boundaryNodes[i];
            

            K.push_back(Trip(i,j,val))
        }
    }

} // namespace std


#endif	