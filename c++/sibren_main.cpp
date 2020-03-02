#include "sibren_functions.hpp"
#include "constants.hpp"
#include <cassert>
#include <cmath>
#include <vector>

namespace pear {



// TODO: fix mesh hier
mesh originalMesh;

// -----------------------------------------------------------------------------


int numberElements; // Number elements in mesh

// Allocate zero-vectors:
std::vector<double> K(pow(2*numberElements,2), 0.0); // Stiffness matrix (will be sparse...)
std::vector<double> K(2*numberElements, 0.0); // Constant term f


void createLinearSystem(originalMesh, K, f);











}
