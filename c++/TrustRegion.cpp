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
Central finite difference approximation of the Jacobian.

INPUT
- Ffun: objective function F(x) of which the Jacobian J(x0) is approximated
- input arguments of objective funnction: std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel
- x0 (column vector) this the position of the current iteration

OUTPUT
- f0 = F(x0) (column vector)
- J = Jacobian
*/
void finite_difference_jacob(arrayxd &f0, spmat & J, void (*Ffun)(spmat &, arrayxd&, arrayxd &, arrayxd&, std::mesh&), arrayxd x0, 
        std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel){

    int Nx = x0.size();
    (*Ffun)(Kstelsel, fstelsel, x0, f0, myMesh); // f0 = F(x0)

    // perform finite difference jacobian evaluation
    double h = 1e-6; // stepsize for first order approximation // original 1e-6
    double per2h = 1/(2*h);
    std::vector<Trip> tripletList; 

    arrayxd fplush(x0.size());
    arrayxd fminh(x0.size());
    arrayxd Jcolj(x0.size());

    for (int j = 0; j < Nx; j++){
        arrayxd xplush = x0; // change x0 to reference!
        arrayxd xminh = x0;
        xplush(j) += h;
        xminh(j) -= h;
        (*Ffun)(Kstelsel, fstelsel, xplush, fplush, myMesh);
        (*Ffun)(Kstelsel, fstelsel, xminh, fminh, myMesh);
        Jcolj = (fplush-fminh)*per2h;

        for ( int i = 0; i < Jcolj.size(); i++){
            assert(!isnan(Jcolj(i)));
            if (Jcolj(i) != 0 && !isnan(Jcolj(i))){
                tripletList.push_back(Trip(i, j, Jcolj(i)));
            }
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
    double f = 0.5*(F.cwiseProduct(F).sum()); 
    return f*1e4;
}





void trustRegion(std::mesh &myMesh, arrayxd & x, void (*Ffun)(spmat&, arrayxd&,arrayxd&, arrayxd&, std::mesh&), arrayxd &x0, spmat &Kstelsel, arrayxd &fstelsel){

    // convergence tolerance
    double grad_tol = 2e-20; 
    int max_iters = 200;

    // regularization parameter
    double lambda =1e-16;// 0.5; No dynamical regularisation
    //double lambdascaling = 0.3;
   // double lambdaMin = 1e-16;

    // trust region parameters
    double MAX_Radius = 25;
    double eita = 0.2; // decide on acceptance of step (between 0 and 1/4) step accepted if rho_k > eita
  
    // F = F(x0)
    arrayxd F(x0.size());
    (*Ffun)(Kstelsel, fstelsel,x0, F, myMesh);

    // a log of the iterations
    vectxd grad_iter(max_iters);

    // loop initialization
    x = x0;
    arrayxd xnext = x0;
    int Nx = x0.size();
    int Nf = F.size();
    spmat J(Nf,Nx);
    spmat A(Nx,Nx);
    double stepRadius = 1;
    bool stepNormLimitReached = false;



    // START ITERATING______________________________________________________________________________________________________________________________________________
    for (int k=1; k <= max_iters; k++){
    
        stepNormLimitReached = false;

        // check for divergence
        assert (x.maxCoeff() < 1e6 && x.minCoeff() > -1e6); // equivalent: x.cwiseAbs().maxCoeff() < 1e6 
        
        // evaluate F(x) and it's jacobian J(x)
        finite_difference_jacob(F, J, Ffun, x, myMesh, Kstelsel, fstelsel);

        if (k == 1){
            // write J to matlab file
            matxd Jdense = matxd(J);
            std::ofstream myFileJac;
            myFileJac.open("../matlab/jacobian.m");
            myFileJac<<"Jdense = [ ";
            for (int i = 0; i<Jdense.rows();i++){
                for (int j = 0; j<Jdense.cols();j++){
                    myFileJac<<Jdense(i,j)<<" ";
                }
            }
            myFileJac<<"];";
            myFileJac.close();
            std::cout<<std::endl;
        }

        //convergence criteria
        vectxd grad = J.transpose()*F.matrix(); // gradient of the scalar objective function f(x) // bring initialisation out of loop!
        double inf_norm_grad = grad.cwiseAbs().maxCoeff();

        // store inf_norm_grad in iteration log
        grad_iter(k) = inf_norm_grad;

        // check for convergence
        if (inf_norm_grad < grad_tol ){ 
            grad_iter = grad_iter.head(k); // discards remaining part of iteration log
            return;
        }


    // COMPUTE STEP PK______________________________________________________________________________________________
        A = J.transpose() * J; 
        for (int i=0; i < A.rows(); i++){
            A.coeffRef(i,i) = A.coeffRef(i,i)*(1+lambda);
        }
        A.makeCompressed();
        //Sparse LU solver: 
            // A*pk = b 
            //  A = J'J + lambda * diag(J'J)    (Fletcher)   // bij normale LM is de Hessian A = J'J + lambda * I                
            //  b = -grad;
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solverA;
        solverA.analyzePattern(A);
        solverA.factorize(A);
        if(solverA.info()!=Eigen::Success) {
            std::cout << "trustRegion: error in Eigen Sparse LU factorization" << std::endl;
        }
        vectxd pk = solverA.solve(-grad); // bring initialisation out of loop!

        // check if pk is smaller then stepRadius
        double stepNorm = pk.norm();
        if (stepNorm > stepRadius){
            pk *= (stepRadius/stepNorm);
            stepNormLimitReached = true;
        }
        xnext = x + pk.array(); // x is in array and pk in vector, this should be fixed

    // COMPUTE TRUSTWORTHINESS______________________________________________________________________________________
        // trustworthiness = rho_k = actual reduction / predicted reduction
        double fx = f(Kstelsel,fstelsel, Ffun, x, myMesh); 
        double fxnext = f(Kstelsel,fstelsel, Ffun, xnext, myMesh); 
        double Ared = fx - fxnext; // actual reduction
        double Pred = -( grad.dot(pk) + 0.5* pk.dot(A*pk) )*1e4; // predicted reduction
        double rhok = Ared/Pred;


    // ADAPT stepRadius FOR NEXT ITERATION___________________________________________________________________________
        if (rhok < 0.25){ // bad model, reduce radius
            stepRadius = 0.25 * stepRadius;
        }
        else if (rhok > 0.75 && stepNormLimitReached){ // good model, inrease radius, but not too much
            stepRadius = min(2*stepRadius, MAX_Radius);
           // std::cout<<"          ...UPDATING stepRadius...   "<< stepRadius << std::endl;
        }

    // DECIDE ON ACCEPTANCE OF STEP___________________________________________________________________________
        if ( rhok > eita){
            x = xnext;
           // std::cout<<"          ...UPDATING X...   " << std::endl;
        }

        // travelled distance from x0
        vectxd difference = (x0-x).matrix(); // bring initialisation out of loop!
        double distance = difference.norm();
        

        // print current log  
       std::cout<< std::endl << "End of iteration: "<< k << "  inf_norm_grad = "<< inf_norm_grad << " f(x) = " << fx <<" norm step = "<< pk.norm() <<  " Total travelled distance = "<< distance << std::endl;
       std::cout<< "        Trust region parameters:    rho = "<< rhok <<  " Ared = "<< Ared << " Pred = "<< Pred << "  Trust region radius = " << stepRadius << std::endl;
    }
    
    std::cout<<" MAX_NB_ITERATIONS exceeded"<< std::endl;
    return;
}



} // end namespace std
