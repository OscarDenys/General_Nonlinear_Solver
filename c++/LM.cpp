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
         std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel){

    // assert that gamma and beta are in a reasonable range
    assert(gamma >= 0 && gamma <=1);
    assert(beta >= 0 && beta <=1);

    int max_nb_steps_counter = 1000; // originally 1000

    // initialize t, evaluate function
    double t = 1;
    double f0 = (*fun)(Kstelsel, fstelsel, Ffun, x0, myMesh);         

    trial_x = x0 + t*pk; 

    while ( (*fun)(Kstelsel,fstelsel, Ffun, trial_x, myMesh)*1e-4 > f0*1e-4 + gamma*t*Jpk ){
        // trial step in x
        t = beta*t;
        trial_x = x0 + t*pk;

        // throw error if line search takes too many iterations
        max_nb_steps_counter -= 1;
        //assert(max_nb_steps_counter > 0);
        if (max_nb_steps_counter == 0){
            std::cout << "line_search: MAX NB STEPS EXCEEDED t = "<< t << std::endl;
            return t;
        }
    }
    return t;
}


/*
Finite difference approximation of the Jacobian.

INPUT
- Ffun: F(x) of which the Jacobian J(x0) is approximated
- x0 (column vector) this the position of the current iteration

OUTPUT
- f0 = F(x0) (column vector)
- J = J(x0)
*/
void finite_difference_jacob(arrayxd &f0, spmat & J, void (*Ffun)(spmat &, arrayxd&, arrayxd &, arrayxd&, std::mesh&), arrayxd x0, 
        std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel){

    //J.setZero();
    int Nx = x0.size();
    (*Ffun)(Kstelsel, fstelsel, x0, f0, myMesh); // f0 = F(x0)

    // perform finite difference jacobian evaluation
    double h = 1e-6; // stepsize for first order approximation // original 1e-6
    std::vector<Trip> tripletList; //triplets.reserve(estimation_of_entries); //--> how many nonzero elements in J?

    arrayxd x = x0;
    arrayxd f(x0.size());
    arrayxd Jcolj(x0.size());

    for (int j = 0; j < Nx; j++){
        x(j) += h;
        (*Ffun)(Kstelsel, fstelsel, x, f, myMesh);
        Jcolj = (f-f0)*(1/h);

        for ( int i = 0; i < Jcolj.size(); i++){
            if (Jcolj(i) != 0 && !isnan(Jcolj(i))){
                tripletList.push_back(Trip(i, j, Jcolj(i)));
            }
            //if(isnan(Jcolj(i))){
            //    std::cout<<"fin_diff_J: element is nan for i ="<< i <<" and j = "<< j << std::endl;
            //}
        }
    }
    J.setFromTriplets(tripletList.begin(), tripletList.end()); // the initial content of J is destroyed.
}


/*
Scalar objective function.

INPUT
- Ffun: F(x)
- x
OUTPUT
- f: (double) f(x) = 0.5*L2-norm(F(x))
*/
double f(spmat & Kstelsel, arrayxd& fstelsel, void (*Ffun)(spmat &, arrayxd &, arrayxd&, arrayxd&, std::mesh&), arrayxd &x, std::mesh &myMesh){
    arrayxd F(x.size());
    (*Ffun)(Kstelsel, fstelsel,x, F, myMesh);
    double f = 0.5*(F.cwiseProduct(F).sum()); // dit is OK
    return f*1e4;
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
void minimize_lm(std::mesh &myMesh, arrayxd & x, void (*Ffun)(spmat&, arrayxd&,arrayxd&, arrayxd&, std::mesh&), arrayxd &x0, spmat &Kstelsel, arrayxd &fstelsel){

    // convergence tolerance
    double grad_tol = 1e-8; // original 1e-4 
    int max_iters = 200;
  
    arrayxd F(x0.size());
    (*Ffun)(Kstelsel, fstelsel,x0, F, myMesh);
    int Nx = x0.size();
    int Nf = F.size();

    // a log of the iterations
    //matxd x_iter(Nx, max_iters);
    vectxd grad_iter(max_iters);

    // line search parameters
    double gamma = 0.01; 
    double beta = 0.8; // original 0.6

    // Levenberg-Marquardt: norm penalisation parameter lambda
    double lambda = 300;//0.5; // lamda large causes the step to be more like gradient, lambda small causes the step to be more like Gauss-Newton
    double lamdascaling = 0.2; 


    // loop initialization
    x = x0;
    spmat J(Nf,Nx);
    spmat A(Nx,Nx);

    //int count_fxsmall = 0;

    for (int k=1; k <= max_iters; k++){

        // check for divergence
        assert (x.maxCoeff() < 1e6 && x.minCoeff() > -1e6); // equivalent: x.cwiseAbs().maxCoeff() < 1e6

        // evaluate F and it's jacobian J
        finite_difference_jacob(F, J, Ffun, x, myMesh, Kstelsel, fstelsel);

        //convergence criteria
        vectxd grad = J.transpose()*F.matrix(); // gradient of the scalar objective function f(x)
        double inf_norm_grad = grad.cwiseAbs().maxCoeff();

        // store x_k and inf_norm_grad in iteration log
        //x_iter.col(k) = x; 
        grad_iter(k) = inf_norm_grad;

        // check for convergence
        if (inf_norm_grad < grad_tol ){
            //x_iter = x_iter.block(0,0,Nx,k); // this gave segmentation fault
            grad_iter = grad_iter.head(k);
            //std::cout<< "solution : "<< x  << std::endl;
            return;
        }

        // find the search direction pk 

        A = J.transpose() * J; 
        for (int i=0; i < A.rows(); i++){
            //std::cout<<"minimize_lm: A.coeffRef(i,i) ="<< A.coeffRef(i,i)<< std::endl;
            A.coeffRef(i,i) = A.coeffRef(i,i)  + lambda;
           // std::cout<<"minimize_lm: A.coeffRef(i,i) ="<< A.coeffRef(i,i)<< std::endl;
        }
        
        A.makeCompressed();
        //Sparse LU solver: 
            // A*pk = b 
            //  A = J'J + lambda * I                //scale invariant version (Fletcher) : J'J + lambda * diag(J'J) ; **zie uitleg onderaan
            //  b = -grad;
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solverA;
        solverA.analyzePattern(A);
        solverA.factorize(A);
        if(solverA.info()!=Eigen::Success) {
            std::cout << "minimize_lm: error in Eigen Sparse LU factorization" << std::endl;
        }
        vectxd pk = solverA.solve(-grad); 
        
        // line search: x = x_new
        double Jpk = grad.dot(pk); 
       // line_search(x, f, Ffun, x, Jpk, pk, gamma, beta, myMesh,  Kstelsel, fstelsel);

        if (k != 1) {
            arrayxd xtrial = x;
            double fxcurrent = f(Kstelsel,fstelsel, Ffun, x, myMesh);
            double t = line_search(xtrial, f, Ffun, x, Jpk, pk, gamma, beta, myMesh,  Kstelsel, fstelsel);
            double fxtrial = f(Kstelsel,fstelsel, Ffun, xtrial, myMesh);
        //if (k > 1) {
            if (fxtrial < fxcurrent){ // && lambda > 0.3){ // step accepted and lambda decreased to make the next step more like the Gauss Newton step (faster than gradient descent)
                x = xtrial;
                if (lambda > 15){
                    lambda = lambda * lamdascaling;
                    }
                std::cout<< "            fxcurrent = " << fxcurrent <<" fxtrial = "<<fxtrial << std::endl;
                std::cout<<"UPDATE X"<< std::endl;
            }
           // if (fxtrial == fxcurrent){
            //    lambda = lambda / lamdascaling;
            //}
            else {//if (fxcurrent <= fxtrial) { // step declined and lambda increased to make the next step more like the stable gradient descent step
                std::cout<< "            fxcurrent = " << fxcurrent <<" fxtrial = "<<fxtrial << std::endl;
                //std::cout<<"LAMBDA INCREASE"<< std::endl;
                lambda = lambda / lamdascaling;
            }

            //if( isnan(fxcurrent)){
            //    std::cout<< "            x current = " << x << std::endl;
            //    std::cout<< "            x trial = " << xtrial << std::endl;
            //    
            //}
            //if (fxtrial < 1e-16) {
            //        std::cout<<"minimize_lm: fxnew smaller than 1e-16 counter ="<< count_fxsmall<< std::endl;
            //        count_fxsmall = count_fxsmall + 1;
            //}
    
        }
        

        // print current log
        std::cout<< "end of iteration: "<< k << "  inf_norm_grad = "<< inf_norm_grad << " f(x) = " << f(Kstelsel,fstelsel, Ffun, x, myMesh) <<" norm step = "<< pk.norm() << " lambda = " << lambda << " x(1) = "<< x(1) << std::endl;

        //if (count_fxsmall ==2){
            //std::cout<< " f(x) : "<< f(Kstelsel,fstelsel, Ffun, x, myMesh)  << std::endl; // DIT GEEFT EEN ERROR
            //std::cout<< " f(x0) : "<< f(Kstelsel,fstelsel, Ffun, x0, myMesh)  << std::endl;
            /* Assertion `row>=0 && row<rows() && col>=0 && col<cols()' failed komt uit Eigen/src/SparseCore/SparseMatrix.h:208: */
        //    return;
        ///}
    }
    
    std::cout<<"minimize_lm: MAX_NB_ITERATIONS exceeded"<< std::endl;
    return;
}


} // end namespace std

/*  **Uitleg scale invariant version (Fletcher) : J'J + lambda * diag(J'J) 

Levenberg's algorithm has the disadvantage that if the value of damping factor lambda  is large, inverting J'J + lambda * I is not used at all.
 Fletcher provided the insight that we can scale each component of the gradient according to the curvature, 
  so that there is larger movement along the directions where the gradient is smaller. 
 This avoids slow convergence in the direction of small gradient. 
 Therefore, Fletcher in his 1971 paper A modified Marquardt subroutine for non-linear least squares replaced the identity matrix
 with the diagonal matrix consisting of the diagonal elements of J'J, thus making the solution scale invariant:*/