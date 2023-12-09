#include "solver.h"
#include "verification.h"
#include "output.h"
// Convenience timer macros for begin/end functions
#define FUNC_BEGIN_TIMER gt.BeginTimer(__func__);
#define FUNC_END_TIMER gt.EndTimer (__func__);
// GRVY Timer object
GRVY::GRVY_Timer_Class gt; // GRVY Timer

solver::solver(int order, int x_mesh, int y_mesh, int solver_type, std::vector<double> q_vector):
    _order(order), xmesh(x_mesh), ymesh(y_mesh), _solver_type(solver_type), q_vec(q_vector){
        // if (order==4){assert(solver_type!=1);}
        assert(order==4 || order==2);

        // xmesh++;
        // if (ymesh!=1){ymesh++;}
        // n_unknowns=xmesh*ymesh;
        if(ymesh!=1) { n_unknowns = (xmesh+1)*(ymesh+1);}
        else { n_unknowns = xmesh + 1;}
        // for now just setting k and q_val on init
        k=1;
        q_val=0;
        // boundary conditons hard coded below - probably worth changing
        
        bc=0;
}
void solver::_1D(){
    FUNC_BEGIN_TIMER
    PetscInt nnz;
    double coef; // coefficient outfront
    // should be done from initialization ?
    if (_order==2){
        // number of nonzeros per row
        nnz=3;
        coef=-k*(xmesh-1)*(xmesh-1); //coefficent out front for 2nd order
    }else if(_order==4){
        nnz=5;
        coef=-k*(xmesh-1)*(xmesh-1)/12;
    }
    // b rhs vector
    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b,PETSC_DECIDE, n_unknowns);
    VecSetFromOptions(b);

    // Matrix A
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_unknowns, n_unknowns, nnz, NULL, &A);
    MatSetFromOptions(A);
    MatSetUp(A);

    // boundary conditions T_b(0)
    MatSetValue(A, 0, 0, 1, INSERT_VALUES); //setting A_00 to 1
    VecSetValue(b, 0, bc,  INSERT_VALUES); //setting b_0 to the boundary condition
    // T_b(n)
    MatSetValue(A, n_unknowns-1,n_unknowns-1, 1, INSERT_VALUES); //setting A_n-1n-1 to 1
    VecSetValue(b, n_unknowns-1, bc,  INSERT_VALUES); //setting b_n-1 to the boundary condition

    // q vector
    std::vector<PetscScalar> q(n_unknowns);
    q=q_vec;
    for (int i=1; i<n_unknowns-1; i++){
        PetscInt idx=i;
        if (nnz==3){ //2nd order
            const PetscInt cols[]={i-1, i, i+1};
            const PetscScalar vals[]={coef, -2*coef, coef};
            MatSetValues(A, 1, &idx, 3 , cols, vals, INSERT_VALUES);
        }else if (nnz==5 && i!=1 && i!=n_unknowns-2){ //4th order
            const PetscInt cols[]={i-2, i-1, i, i+1, i+2};
            const PetscScalar vals[]={-1*coef, 16*coef, -30*coef, 16*coef, -1*coef};
            MatSetValues(A, 1, &idx, nnz , cols, vals, INSERT_VALUES);
        }
        else if (nnz == 5 && i == 1){
            const PetscInt cols[]={1, 2, 3, 4, 5};
            const PetscScalar vals[]={12*coef, -48*coef, 72*coef, -48*coef, 12*coef};
            MatSetValues(A, 1, &idx, nnz , cols, vals, INSERT_VALUES);
        }
        else if (nnz ==5 && i == n_unknowns - 2){
            const PetscInt cols[]={i-4, i-3, i-2, i-1, i};
            const PetscScalar vals[]={12*coef, -48*coef, 72*coef, -48*coef, 12*coef};
            MatSetValues(A, 1, &idx, nnz , cols, vals, INSERT_VALUES);
        }
        //can be set with all of q vetcor at once instead of in for loop - might save time but idk
        VecSetValues(b, 1, &idx, &q[i],  INSERT_VALUES);
    }
    //SOLVE
    _setup_solve_destroy(A, b);
    FUNC_END_TIMER
    //gt.Summarize();
}

void solver::_2D(){
    FUNC_BEGIN_TIMER
    PetscInt nnz; // number of nonzeros per row
    // should maybe be done from initialization?
    if (_order==2){
        nnz=5;
    }else if(_order==4){
        nnz=9;
    }
    // b rhs vector
    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b,PETSC_DECIDE, n_unknowns);
    VecSetFromOptions(b);

    // Matrix A
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_unknowns, n_unknowns, nnz, NULL, &A);
    MatSetFromOptions(A);
    MatSetUp(A);
    int total=0;
    std::vector<PetscScalar> q(n_unknowns);
    q=q_vec;

    // BOUNDARY CONDITIONS
    for (int i=0; i<xmesh+1; i++){
        PetscInt idx=i;
        //setting T(x,y=0)=bc_0
        MatSetValue(A, idx, idx, 1, INSERT_VALUES); //setting first xmesh diagonals of A to 1
        VecSetValue(b, idx, bc,  INSERT_VALUES); //setting bc_0 to the first xmesh entries of b
        //setting T(x,y=1)=bc_1
        MatSetValue(A, n_unknowns-idx-1,n_unknowns-idx-1, 1, INSERT_VALUES); //setting last xmesh diagonals of A to 1
        VecSetValue(b, n_unknowns-1-idx, bc,  INSERT_VALUES); //setting bc_1 to the first xmesh entries of b
        //setting T(x=0,y) = T(x=1,y)=0
        if(i!=0 && i!=xmesh){
            idx=i*(ymesh+1);
            MatSetValue(A, idx,idx, 1, INSERT_VALUES); //x=0
            MatSetValue(A, idx+ymesh,idx+ymesh, 1, INSERT_VALUES);//x=1
            VecSetValue(b, idx, bc, INSERT_VALUES);
            VecSetValue(b, idx+ymesh,bc, INSERT_VALUES);
        }
    }
    double x_2=(xmesh)*(xmesh);
    double y_2=(ymesh)*(ymesh);
    for (int i=xmesh+1; i<n_unknowns-xmesh-1; i++){
        if (i%(ymesh+1)!=0 && (i+1)%(ymesh+1)!=0){ // this if statement forces loop not to overwrite x=0,1 boundary conditions
            PetscInt idx=i;
            if (nnz==5 || i<=2*(xmesh+1) || i>=n_unknowns-2*(xmesh+1) || i%(ymesh+1)==1 ||i%(ymesh+1)==ymesh-1 ){ //2nd order
                // 4th uses 2nd order for coeffiecents next to boundary conditions
                const PetscInt cols[]={i-xmesh-1,i-1, i, i+1, i+xmesh+1};
                const PetscScalar vals[]={-k*y_2, -k*x_2, 2*k*(x_2+y_2), -k*x_2, -k*y_2};
                MatSetValues(A, 1, &idx, 5 , cols, vals, INSERT_VALUES);
            }else if (nnz==9){ //4th order
                const PetscInt cols[]={i-2*(xmesh+1), i-xmesh-1, i-2, i-1, i, i+1, i+2, i+xmesh+1, i+2*(xmesh+1)};
                const PetscScalar vals[]={k*y_2/12, -16*k*y_2/12, k*x_2/12, -16*k*x_2/12,
                            30*k*(x_2+y_2)/12, -16*k*x_2/12, k*x_2/12, -16*k*y_2/12, k*y_2/12};
                MatSetValues(A, 1, &idx, nnz , cols, vals, INSERT_VALUES);
            }
            //can be set with all of q vetcor at once instead of in for loop - might save time but idk
            VecSetValues(b, 1, &idx, &q[i],  INSERT_VALUES);
        }
    }
    //SOLVE
    _setup_solve_destroy(A, b);
    FUNC_END_TIMER
    //gt.Summarize();
}
void solver::_setup_solve_destroy(Mat A, Vec b){
    FUNC_BEGIN_TIMER
    // solution vector x
    VecDuplicate(b, &x);

    //finish set up on matrices
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);


    // Displaying the matrices
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "output_A.txt", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "output_B.txt", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
    VecView(b, viewer);
    PetscViewerDestroy(&viewer);


    //creating and setting up the solve
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetTolerances(ksp, 1e-20, 1e-20 , PETSC_DEFAULT  , _maxit);

    // will use GMRES solve with ILU preconditioning (petsc default ksp solve)
    // if solver type not set to 1 or 2

    if(_solver_type==1){ // jacobi
        KSPSetType(ksp, KSPRICHARDSON);
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCJACOBI);
    }else if(_solver_type==2){ //gauss siedel
        KSPSetType(ksp, KSPPREONLY);
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCSOR);
        PCSORSetIterations(pc, 100, 10);
    }
    //solve
    KSPSolve(ksp, b, x);

    PetscInt size;
    VecGetSize(x, &size);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "output_x.txt", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);

    //getting solution data in c++ vector
    std::vector<double> solution_vector(size);
    PetscScalar* petsc_data;
    VecGetArray(x, &petsc_data);
    std::copy(petsc_data, petsc_data + size, solution_vector.begin());

    KSPDestroy(&ksp);
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);

    populate_file(xmesh, ymesh, solution_vector);
    _temperature = solution_vector;
    FUNC_END_TIMER
    //gt.Summarize();
}

std::vector<double> solver::solve(double tol, int maxit){
    _tol=tol;
    _maxit=maxit;
    if (ymesh==1){
        _1D();
    }else{
        assert(ymesh>1);
        _2D();
    }
    FUNC_END_TIMER
    gt.Summarize();
    return _temperature;
}
