#ifndef SOLVER_H
#define SOLVER_H

#include <petsc.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include<grvy.h>
#include<sys/time.h>
#include<time.h>
#include<unistd.h>
#include <cmath>

class solver{
    public:
        /**
        *@brief default constructor for solver class.
        **/
        solver(){}

        /**
        *@brief constructor for solver class.
        *@param order: 2nd or 4th order - should either equal 2 or 4 exactly
        *@param x_mesh: discretations in x direction
        *@param y_mesh: discretations in y direction NOTE if y_mesh ==1 1D solve inferred
        *@param solver_type: int -- 1 if Jacobi -- 2 if gauss-seidel -- 3 for gmres
        *@param q_vector : q vector aka rhs vector
        **/
        solver(int order, int x_mesh, int y_mesh, int solver_type, std::vector<double> q_vector);

        /**
        *@brief solves linear system according to constructor inputs
        **/
        std::vector<double> solve(double tol=1e-12, int maxit=1e5);



    private:
        /**
        *@brief either 2nd(_order==2) or 4th(_order==4) order
        **/
        int _order;

        /**
        *@brief number of discretations in x direction
        **/
        int xmesh;

        /**
        *@brief number of discretations in y direction
        **/
        int ymesh;

        /**
        *@brief -- 1 if Jacobi -- 2 if gauss-seidel -- 3 for gmres
        **/
        int _solver_type;

        /**
        *@brief  a vector of temperatures [K] (solution vector)
        **/
        std::vector<double> _temperature;

        /**
        *@brief Matrix A
        **/
        Mat A;

        /**
        *@brief  b is the rhs vector and x is the solution vector
        **/
        Vec b,x;

        /**
        *@brief krylov subspace method object
        **/
        KSP ksp;

        /**
        *@brief preconditioner object
        **/
        PC pc;

        /**
        *@brief the 1D solver
        **/
        void _1D();

        /**
        *@brief the 2D solver
        **/
        void _2D();

        /**
        *@brief Finalizes matrix preparation, GMRES solve with ILU preconditioning (petsc default ksp solve), and then destroys matrix A, vector b, and solution vector x
        *@param A lhs PetscMatrix A
        *@param b rhs Petsc vector b
        **/
        void _setup_solve_destroy(Mat A, Vec b);

        /**
        *@brief maximum number of iterations
        **/
        int _maxit;

        /**
        *@brief relative tolerance for convergence
        **/
        int _tol;

        /**
        *@brief thermal conductivity [W/(m*K)]
        **/
        double k;

        /**
        *@brief number of unknowns in system
        **/
        PetscInt n_unknowns;

        /**
        *@brief maybe just for testing but currently setting rhs to q_val at every entry other than boundary conditions
        **/
        double q_val;

        /**
        *@brief boundary condition at all edges = 0
        **/
        PetscScalar bc;

        /**
        *@brief q vector to be passed in
        **/
        std::vector<double> q_vec;
};
#endif