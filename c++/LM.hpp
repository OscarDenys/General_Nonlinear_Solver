#include <cassert>
#include <cmath>
#include <vector>
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include <iostream>
#include "mesh.hpp"


typedef Eigen::SparseMatrix<double> spmat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;
typedef Eigen::VectorXd vectxd;
typedef Eigen::MatrixXd matxd;
typedef Eigen::ArrayXd arrayxd;


namespace std {

/*
Line search using Armijo conditions and backtracking.

INPUT
- fun: scalar function f = 0.5*L2-norm(F(x))
- F: vectxd F = F(x) (input for f(F) = 0.5*L2-norm(F))
- initial guess x0
- Jpk: J is gradient(!) of fun at x0, Jpk is J'*pk (scalar)
- pk: search direction
- gamma: armijo condition scaling function
- beta: backtracking parameter

OUTPUT
- trial_x = x0 + t*pk
- t: scaling of the step pk (returned)
*/
// TODO misschien argument x0 en x samen gebruiken
double line_search(arrayxd & trial_x, double (*fun)(spmat&, arrayxd&, void (*Ffun)(spmat &, arrayxd &, arrayxd&, arrayxd&, std::mesh&), arrayxd& x, 
        std::mesh&), void (*Ffun)(spmat &, arrayxd &, arrayxd&, arrayxd&, std::mesh&),  arrayxd& x0, double Jpk, arrayxd pk, double gamma, double beta,
         std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel);


/*
Finite difference approximation of the Jacobian.

INPUT
- Ffun: F(x) of which the Jacobian J(x0) is approximated
- x0 (column vector)

OUTPUT
- f0 = F(x0) (column vector)
- J = J(x0)
*/
void finite_difference_jacob(arrayxd &f0, spmat & J, void (*Ffun)(spmat &, arrayxd&, arrayxd &, arrayxd&, std::mesh&), arrayxd x0, 
        std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel);

/*
Scalar objective function.

INPUT
- Ffun: F(x)
- x
OUTPUT
- f: (double) f(x) = 0.5*L2-norm(F(x))
*/
double f(spmat & Kstelsel, arrayxd& fstelsel, void (*Ffun)(spmat &, arrayxd &, arrayxd&, arrayxd&, std::mesh&), arrayxd &x, std::mesh &myMesh);


/*
Levenberg-Marquardt algorithm.

INPUT
- Ffun: function F(x) to be minimized
- starting point vector x0
- alpha_k: weighting parameter to penalize norm(step) (TODO)
OUTPUT
- solution vector x
- x_iter: each of the intermediate values xk
- grad_iter: norm of gradient in each iteration
*/
// TODO (eventueel) :void minimize_lm(vectxd x, matxd x_iter, vectxd grad_iter, arrayxxd (*Ffun)(vectxd), vectxd x0);
void minimize_lm(std::mesh &myMesh, arrayxd & x, void (*Ffun)(spmat&, arrayxd&,arrayxd&, arrayxd&, std::mesh&), arrayxd &x0, spmat &Kstelsel, arrayxd &fstelsel);


} // end namespace std