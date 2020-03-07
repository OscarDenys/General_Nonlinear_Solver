#ifndef helpers_hpp
#define helpers_hpp
#include <vector>
#include "Eigen/SparseCore"

namespace std {

    double detJac(vector<float> P1, vector<float> P2, vector<float> P3);
    void applySigmaAndAddCommonPart(std::vector<double> &result, std::vector<double> const &commonPart);

    void evaluateRespiration(int nodeIndex, Eigen::VectorXd &prevSol, double &Ru, double &Rv);
    void evaluateRespiration(int nodeIndex1, int nodeIndex2, Eigen::VectorXd &prevSol, double &Ru, double &Rv);
    double evaluateRu(double Cu, double Cv);
    double evaluateRv(double Cu, double Cv, double Ru);

} // namespace std


#endif