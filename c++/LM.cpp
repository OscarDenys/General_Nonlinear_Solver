#include <cassert>
#include <cmath>
#include <vector>
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include <iostream>

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
double line_search(arrayxd & trial_x, double (*fun)(arrayxd), arrayxd F,  arrayxd x0, double Jpk, arrayxd pk, double gamma, double beta){

    // assert that gamma and beta are in a reasonable range
    assert(gamma >= 0 && gamma <=1);
    assert(beta >= 0 && beta <=1);

    int max_nb_steps_counter = 1000;

    // initialize t, evaluate function
    double t = 1;
    double f0 = (*fun)(F);

    trial_x = x0 + t*pk;
    while ( (*fun)(trial_x) > f0 + gamma*t*Jpk ){
        // trial step in x
        t = beta*t;
        trial_x = x0 + t*pk;

        //std::cout << "(*fun)(trial_x) = " << (*fun)(trial_x) << std::endl; 
        //std::cout << "f0 + gamma*t*Jpk = " << f0 + gamma*t*Jpk << std::endl;

        //std::cout << "trial_x = " << trial_x << std::endl; // NIET OK -- iets met Jpk?


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
void finite_difference_jacob(arrayxd & f0, spmat & J, arrayxd (*Ffun)(arrayxd), arrayxd x0){
/*
    // make sure x0 is a column vector 
    if (x0.cols() == 1){
        std::cout<< "finite_difference_jacob: x0 needs to be column vector"<< std::endl;
        x0.transpose(); 
    }
*/
    int Nx = x0.size();
    f0 = (*Ffun)(x0);
/*
    // make sure fun returns a column vector
    f0 = (*Ffun)(x0);
    if (f0.cols() == 1){
        std::cout<< "finite_difference_jacob: fun needs to return a column vector"<< std::endl;
        f0.transpose();
    }
    //int Nf = f0.rows();

    // initialize empty J (this would be bad practice memory-wise..)
    //J = spmat(Nf,Nx);
*/
    // perform finite difference jacobian evaluation
    double h = 1e-6; // stepsize for first order approximation
    std::vector<Trip> tripletList; //triplets.reserve(estimation_of_entries); //--> how many nonzero elements in J?

    for (int j = 0; j < Nx; j++){
        arrayxd x = x0;
        x(j) += h;
        arrayxd f = (*Ffun)(x);

        //std::cout << "f = " << f << endl; OK

        arrayxd Jcolj = (f-f0)*(1/h);

        for ( int i = 0; i < Jcolj.size(); i++){
            if (Jcolj(i) != 0){
                tripletList.push_back(Trip(i, j, Jcolj(i)));
            }
        }
    }
    J.setFromTriplets(tripletList.begin(), tripletList.end());

    //std::cout << "J = " << matxd(J) << std::endl; OK
}

/*
Scalar objective function.

INPUT
- Ffun: F(x)
- x
OUTPUT
- f: (double) f(x) = 0.5*L2-norm(F(x))
*/
/*
double f(vectxd (*Ffun)(vectxd), vectxd x){
    vectxd F = (*Ffun)(x);
    double f = 0.5*(F.dot(F));
    return f;
}
*/

/*
Scalar objective function.

INPUT
- F = F(x) uitgewerkt voor een zekere x
OUTPUT
- f: (double) f(x) = 0.5*L2-norm(F)
*/
double f(arrayxd F){
    double f = 0.5*(F.cwiseProduct(F).sum());
    //std::cout << "f = " << f << std::endl;  // OK
    //std::cout << "F = " << F << std::endl;
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
//void minimize_lm(vectxd x, matxd x_iter, vectxd grad_iter, arrayxxd (*Ffun)(vectxd), vectxd x0);
void minimize_lm(arrayxd & x, arrayxd (*Ffun)(arrayxd), arrayxd x0){

    // convergence tolerance
    double grad_tol = 1e-4;
    int max_iters = 200;
/*
    // make sure x0 is a column vector 
    if (x0.cols() == 1){
        std::cout<< "minimize_lm: x0 needs to be column vector"<< std::endl;
        x0.transpose(); 
    }
    int Nx = x0.rows();

    // make sure fun returns a column vector
    arrayxd F = (*Ffun)(x0);
    if (F.cols() == 1){
        std::cout<< "minimize_lm: fun needs to return a column vector"<< std::endl;
        F.transpose();
    }
    int Nf = F.rows();
*/    
    arrayxd F = (*Ffun)(x0);
    int Nx = x0.size();
    int Nf = F.size();

    // a log of the iterations
    matxd x_iter(Nx, max_iters);
    vectxd grad_iter(max_iters);

    // line search parameters
    double gamma = 0.01; 
    double beta = 0.6;

    // loop initialization
    x = x0;
    spmat J(Nf,Nx);

    for (int k=1; k <= max_iters; k++){

        // check for divergence
        assert (x.maxCoeff() < 1e6 && x.minCoeff() > -1e6); // or x.cwiseAbs().maxCoeff() < 1e6

        // evaluate F and it's jacobian J
        finite_difference_jacob(F, J, Ffun, x);

        //convergence criteria

        //std::cout << "F = " << F << std::endl; // OK
        //std::cout << "J matrix = " << matxd(J) << std::endl; // steeds de eerste J.. ---> steeds dezelfde x --> while conditie in line search

        vectxd grad = J.transpose()*F.matrix(); // gradient of the scalar objective function f(x)

        //std::cout << "grad = " << grad << std::endl; // NIET OK

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

        // find the search direction
        spmat A = J.transpose() * J;
        A.makeCompressed();
        //Sparse LU solver (square system pk = -(J'*J)\(J'*F) )
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solverA;
        solverA.analyzePattern(A);
        solverA.factorize(A);
        if(solverA.info()!=Eigen::Success) {
            std::cout << "minimize_lm: error in Eigen Sparse LU factorization" <<"\n";
        }
        vectxd pk = solverA.solve(-grad); // solver is OK

        std::cout << "pk = " << pk << std::endl; // 2de pk is fout

        //TODO: A+ Identity_matrix * weight parameter alpha_k! // matrix inversion!!! (this is the Gauss-Newton method without alpha_k)
        // b = -grad; A = J'J 
        // NOTE: A should be in compressed and column major order form! (A*pk=b)
        // SOURCES for code sparse solver:
        // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
        // https://scicomp.stackexchange.com/questions/21343/solving-linear-equations-using-eigen
        
        // line search
        double Jpk = grad.dot(pk); 
        line_search(x, f, F, x, Jpk, pk, gamma, beta);
        //std::cout << "x = " << x << std::endl; // NOT CHANGING
    }
    
    std::cout<<"minimize_lm: MAX_NB_ITERATIONS exceeded"<< std::endl;
    return;
}


} // end namespace std