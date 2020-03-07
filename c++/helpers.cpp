#include <vector>
#include <cmath>
#include <cassert>
#include "helpers.hpp"
#include "constants.hpp"
#include "Eigen/SparseCore"


namespace std {

    double detJac(vector<float> P1, vector<float> P2, vector<float> P3) {
        return double(P1[1]*P2[0] - P1[0]*P2[1] + P1[0]*P3[1] - P1[1]*P3[0] - P2[0]*P3[1] + P2[1]*P3[0]); // Bewerking op single prec floats opslaan in double precision float mag???
    } // detJac

    void applySigmaAndAddCommonPart(std::vector<double> &result, std::vector<double> const &commonPart){
        // TODO: volgens wiskundige afleiding moet sigma_uz maal r gedeelte en andersom, dit is niet hoe wij het in MATLAB doen!!!
        result[0] = sigma_uz*commonPart[0] + sigma_ur*commonPart[1];
        result[1] = sigma_vz*commonPart[0] + sigma_vr*commonPart[1];

        return;
    } // applySigmaAndAddCommonPart()

    // respiration kinetics: R - formula (3)
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

} // namespace std