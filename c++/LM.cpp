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
double line_search(arrayxd & trial_x, double (*fun)(void (*Ffun)(arrayxd, arrayxd, std::mesh&), arrayxd x, 
        std::mesh&), void (*Ffun)(arrayxd, arrayxd, std::mesh&),  arrayxd x0, double Jpk, arrayxd pk, double gamma, double beta, std::mesh &myMesh){

    // assert that gamma and beta are in a reasonable range
    assert(gamma >= 0 && gamma <=1);
    assert(beta >= 0 && beta <=1);

    int max_nb_steps_counter = 1000;

    // initialize t, evaluate function
    double t = 1;
    double f0 = (*fun)(Ffun, x0, myMesh);         

    trial_x = x0 + t*pk; 

    while ( (*fun)(Ffun, trial_x, myMesh) > f0 + gamma*t*Jpk ){
        // trial step in x
        t = beta*t;
        trial_x = x0 + t*pk;

        // throw error if line search takes too many iterations
        max_nb_steps_counter -= 1;
        assert(max_nb_steps_counter > 0);
    }

    return t;
}


/*
Finite difference approximation of the Jacobian.

INPUT
- Ffun: F(x) of which the Jacobian J(x0) is approximated
- x0 (column vector)

OUTPUT
- f0 = F(x0) (column vector)
- J = J(x0)
*/
void finite_difference_jacob(arrayxd &f0, spmat & J, void (*Ffun)(arrayxd, arrayxd, std::mesh&), arrayxd x0, std::mesh &myMesh){

    int Nx = x0.size();
    (*Ffun)(x0, f0, myMesh);

    // perform finite difference jacobian evaluation
    double h = 1e-6; // stepsize for first order approximation
    std::vector<Trip> tripletList; //triplets.reserve(estimation_of_entries); //--> how many nonzero elements in J?

    arrayxd x = x0;
    arrayxd f(x0.size());
    arrayxd Jcolj(x0.size());

    for (int j = 0; j < Nx; j++){
        x(j) += h;
        (*Ffun)(x, f, myMesh);
        Jcolj = (f-f0)*(1/h);

        for ( int i = 0; i < Jcolj.size(); i++){
            if (Jcolj(i) != 0){
                tripletList.push_back(Trip(i, j, Jcolj(i)));
            }
        }
    }
    J.setFromTriplets(tripletList.begin(), tripletList.end());
}


/*
Scalar objective function.

INPUT
- Ffun: F(x)
- x
OUTPUT
- f: (double) f(x) = 0.5*L2-norm(F(x))
*/
double f(void (*Ffun)(arrayxd, arrayxd, std::mesh&), arrayxd x, std::mesh &myMesh){
    arrayxd F(x.size());
    (*Ffun)(x, F, myMesh);
    double f = 0.5*(F.cwiseProduct(F).sum()); // dit is OK
    return f;
}


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
void minimize_lm(std::mesh &myMesh, arrayxd & x, void (*Ffun)(arrayxd, arrayxd, std::mesh&), arrayxd x0){

    // convergence tolerance
    double grad_tol = 1e-4;
    int max_iters = 200;
  
    arrayxd F(x0.size());
    (*Ffun)(x0, F, myMesh);
    int Nx = x0.size();
    int Nf = F.size();

    // a log of the iterations
    matxd x_iter(Nx, max_iters);
    vectxd grad_iter(max_iters);

    // line search parameters
    double gamma = 0.01; 
    double beta = 0.6;

    // Levenberg-Marquardt: norm penalisation parameter lambda
    double lambda = 0.01;

    // loop initialization
    x = x0;
    spmat J(Nf,Nx);

    for (int k=1; k <= max_iters; k++){

        // check for divergence
        assert (x.maxCoeff() < 1e6 && x.minCoeff() > -1e6); // or x.cwiseAbs().maxCoeff() < 1e6

        // evaluate F and it's jacobian J
        finite_difference_jacob(F, J, Ffun, x, myMesh);

        //convergence criteria
        vectxd grad = J.transpose()*F.matrix(); // gradient of the scalar objective function f(x)
        double inf_norm_grad = grad.cwiseAbs().maxCoeff();

        // store x_k and inf_norm_grad in iteration log
        x_iter.col(k) = x;
        grad_iter(k) = inf_norm_grad;

        // print current log
        std::cout<< "iteration: "<< k << "  inf_norm_grad = "<< inf_norm_grad << std::endl;

        // check for convergence
        if (inf_norm_grad < grad_tol){
            x_iter = x_iter.block(0,0,Nx,k);
            grad_iter = grad_iter.head(k);
            return;
        }

        // find the search direction pk 

        spmat A = J.transpose() * J; 
        for (int i; i < A.rows(); i++){
            A.coeffRef(i,i) = (A.coeffRef(i,i) + 0.01) * lambda;
        }
        
        A.makeCompressed();
        //Sparse LU solver: 
            // A*pk = b 
            //  A = J'J + lambda * (diag(J'J) + 0.01);
            //  b = -grad;
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solverA;
        solverA.analyzePattern(A);
        solverA.factorize(A);
        if(solverA.info()!=Eigen::Success) {
            std::cout << "minimize_lm: error in Eigen Sparse LU factorization" <<"\n";
        }
        vectxd pk = solverA.solve(-grad); 
        
        // line search
        double Jpk = grad.dot(pk); 
        line_search(x, f, Ffun, x, Jpk, pk, gamma, beta, myMesh);
    }
    
    std::cout<<"minimize_lm: MAX_NB_ITERATIONS exceeded"<< std::endl;
    return;
}


} // end namespace std