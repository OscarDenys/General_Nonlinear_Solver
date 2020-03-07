#ifndef sibren_functions_hpp
#define sibren_functions_hpp
#include <cassert>
#include <cmath>
#include <vector>

namespace std {
  void integral1(mesh const &myMesh, std::vector<Eigen::Triplets<double>> &K);
  void applySigmaAndAddCommonPart(std::vector<double> &result, std::vector<double> const &commonPart);
  void integral2nonlinear(mesh & myMesh, Eigen::VectorXd &H, Eigen::VectorXd &C);

}









#endif
