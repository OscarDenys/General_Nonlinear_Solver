#include <cassert>
#include <cmath>
#include <vector>
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include <iostream>
#include "mesh.hpp"

// TODO
// update make file
// make hpp file (check first if first draft is finished)
// test and debug ---> print output after iteration


typedef Eigen::SparseMatrix<double> spmat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;
typedef Eigen::VectorXd vectxd;
typedef Eigen::MatrixXd matxd;
typedef Eigen::ArrayXd arrayxd;


namespace std {


/*
Finite difference approximation of the Jacobian.

INPUT
- Ffun: F(x) of which the Jacobian J(x0) is approximated
- x0 (column vector) this the position of the current iteration

OUTPUT
- f0 = F(x0) (column vector)
- J = J(x0)
*/
/*
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
    int nancount = 0;

    for (int j = 0; j < Nx; j++){
        x(j) += h;
        (*Ffun)(Kstelsel, fstelsel, x, f, myMesh);
        Jcolj = (f-f0)*(1/h);

        for ( int i = 0; i < Jcolj.size(); i++){
            if (Jcolj(i) != 0 && !isnan(Jcolj(i))){
                tripletList.push_back(Trip(i, j, Jcolj(i)));
            }
            if(isnan(Jcolj(i))){
                nancount +=1;
            }
        }
    }
    if (nancount > 0){
        std::cout<<"fin_diff_J: number of nans in jacobian = "<<nancount << std::endl;
    }
    J.setFromTriplets(tripletList.begin(), tripletList.end()); // the initial content of J is destroyed.
}*/

// CENTRAL FINITE DIFFERENCE
void finite_difference_jacob(arrayxd &f0, spmat & J, void (*Ffun)(spmat &, arrayxd&, arrayxd &, arrayxd&, std::mesh&), arrayxd x0, 
        std::mesh &myMesh, spmat &Kstelsel, arrayxd &fstelsel){

    //J.setZero();
    int Nx = x0.size();
    (*Ffun)(Kstelsel, fstelsel, x0, f0, myMesh); // f0 = F(x0)

    // perform finite difference jacobian evaluation
    double h = 1e-6; // stepsize for first order approximation // original 1e-6
    double per2h = 1/(2*h);
    std::vector<Trip> tripletList; //triplets.reserve(estimation_of_entries); //--> how many nonzero elements in J?

    arrayxd fplush(x0.size());
    arrayxd fminh(x0.size());
    arrayxd Jcolj(x0.size());

    for (int j = 0; j < Nx; j++){
        arrayxd xplush = x0; //x0 argument zou in principe reference (&) mogen zijn
        arrayxd xminh = x0;
        xplush(j) += h;
        xminh(j) -= h;
        (*Ffun)(Kstelsel, fstelsel, xplush, fplush, myMesh);
        (*Ffun)(Kstelsel, fstelsel, xminh, fminh, myMesh);
        Jcolj = (fplush-fminh)*per2h;

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





void trustRegion(std::mesh &myMesh, arrayxd & x, void (*Ffun)(spmat&, arrayxd&,arrayxd&, arrayxd&, std::mesh&), arrayxd &x0, spmat &Kstelsel, arrayxd &fstelsel){

    // convergence tolerance
    double grad_tol = 1e-17; // original 1e-4 
    int max_iters = 200;

    // regularization parameter
    double lambda = 0.5;
    double lambdascaling = 0.3;
    double lambdaMin = 1e-16;

    // trust region parameters
    double MAX_Radius = 25;
    double eita = 0.2; // decide on acceptance of step (between 0 and 1/4) step accepted if rho_k > eita
  
    // F = F(x0)
    arrayxd F(x0.size());
    (*Ffun)(Kstelsel, fstelsel,x0, F, myMesh);

    // a log of the iterations
    //matxd x_iter(Nx, max_iters);
    vectxd grad_iter(max_iters);

    // loop initialization
    x = x0;
    arrayxd xnext = x0;
    int Nx = x0.size();
    int Nf = F.size();
    spmat J(Nf,Nx);
    spmat A(Nx,Nx);
    double stepRadius = 1;


    //test objective function
    arrayxd xtest(x0.size());
    std::cout << "f(x0)= " << f(Kstelsel,fstelsel, Ffun, x0, myMesh) << std::endl;
    std::cout << "f(x = 0 )= " << f(Kstelsel,fstelsel, Ffun, xtest, myMesh) << std::endl;
    xtest = 1;
    std::cout << "f(x = 1 )= " << f(Kstelsel,fstelsel, Ffun, xtest, myMesh) << std::endl;
    double nanElem = 0;
    bool stepNormLimitReached = false;



    // START ITERATING______________________________________________________________________________________________________________________________________________
    for (int k=1; k <= max_iters; k++){

        stepNormLimitReached = false;

        // check for divergence
        assert (x.maxCoeff() < 1e6 && x.minCoeff() > -1e6); // equivalent: x.cwiseAbs().maxCoeff() < 1e6 // dit geeft assertion error 
        //if (x.maxCoeff() < 1e6 && x.minCoeff() > -1e6){                                                         // conditional moet nog negatief gezet worden
       //     std::cout << "divergence x "<< x.cwiseAbs().maxCoeff()<< std::endl;
         //   return;
        //}

        // evaluate F(x) and it's jacobian J(x)
        finite_difference_jacob(F, J, Ffun, x, myMesh, Kstelsel, fstelsel);

        //convergence criteria
        vectxd grad = J.transpose()*F.matrix(); // gradient of the scalar objective function f(x) // bring initialisation out of loop
        double inf_norm_grad = grad.cwiseAbs().maxCoeff();

        // store inf_norm_grad in iteration log
        grad_iter(k) = inf_norm_grad;

        // check for convergence
        if (inf_norm_grad < grad_tol ){ 
            grad_iter = grad_iter.head(k); // discard remaining part of iteration log
            return;
        }


    // COMPUTE STEP PK______________________________________________________________________________________________
        A = J.transpose() * J; 
        for (int i=0; i < A.rows(); i++){
            A.coeffRef(i,i) = A.coeffRef(i,i)*(1+lambda);
            // std::cout<<"minimize_lm: A.coeffRef(i,i) ="<< A.coeffRef(i,i)<< std::endl;
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
        vectxd pk = solverA.solve(-grad); // bring initialisation out of loop

        // check if pk is smaller then stepRadius
        double stepNorm = pk.norm();
        if (stepNorm > stepRadius){
            pk *= (stepRadius/stepNorm);
            stepNormLimitReached = true;
            //std::cout<<" BEEP BOOP BIEP BOOP ...stepNorm > stepRadius ...   " << stepNormLimitReached << std::endl;
        }
        // x is in array and pk in vector, this should be fixed
        xnext = x + pk.array(); 

    // COMPUTE TRUSTWORTHINESS______________________________________________________________________________________
        // trustworthiness = rho_k = actual reduction / predicted reduction
        double fx = f(Kstelsel,fstelsel, Ffun, x, myMesh); //misschien helpt dit tegen nan
        double fxnext = f(Kstelsel,fstelsel, Ffun, xnext, myMesh); 
        double Ared = fx - fxnext; // actual reduction
        double Pred = -( grad.dot(pk) + 0.5* pk.dot(A*pk) )*1e4; // predicted reduction
        double rhok = Ared/Pred;

        if (isnan(fx)){
            nanElem = 0;
            for (int i = 0; i < x.size(); i++){
                if ( isnan(x(i)) ){
                    nanElem += 1;
                }
            }
             std::cout << "              f(x) is nan..             number of nan elements in x = " << nanElem << std::endl;
                    //if (nanElem > 0){
                    //    std::cout << "       f(x) is nan..             number of nan elements in x = " << nanElem << std::endl;
                    //}
                    //else{
                    //    std::cout << "       f(x) is nan..              zero nan elements in x" << std::endl;
                    //}
        }

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

        if ( abs(rhok-1) < 0.2 ){
                lambda = min(lambda *lambdascaling, lambdaMin); 
        }
        else if (abs(rhok-1) > 10){
                lambda = lambda/ lambdascaling; 
        }

        // travelled distance from xo
        vectxd difference = (x0-x).matrix(); // bring initialisation out of loop
        double distance = difference.norm();
        

        // print current log  // " Ared = "<< Ared << " Pred = "<< Pred <<
       //std::cout<< "end of iteration: "<< k << "  inf_norm_grad = "<< inf_norm_grad << " f(x) = " << fx <<" norm step = "<< pk.norm()<< " rho = "<< rhok <<  " travelled distance = "<< distance << "  Trust region radius = " << stepRadius << std::endl;
       std::cout<< "        rho = "<< rhok <<  " Ared = "<< Ared << " Pred = "<< Pred << "  Trust region radius = " << stepRadius<< " lambda = "<< lambda << std::endl;
       std::cout<< "end of iteration: "<< k << "  inf_norm_grad = "<< inf_norm_grad << " f(x) = " << fx <<" norm step = "<< pk.norm() <<  " travelled distance = "<< distance << std::endl;
    }
    
    std::cout<<" MAX_NB_ITERATIONS exceeded"<< std::endl;
    return;
}



} // end namespace std
