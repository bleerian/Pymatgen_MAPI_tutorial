#include "solver.h"
#include "input_parse.h"
#include "unit_test.h"


void make_test(){
    for (int solver_type =1; solver_type<4; solver_type++){
        for (int j=1; j<3; j++){ //order
            if (solver_type==1 && j>1){
                continue;
            }
            int order=2*j;
            for (int i =2; i<7; i++){ //mesh size
                int mesh=std::pow(2,i);
                std::cout<<" 1d  " <<order<<"  order  "<<"solver_type=="<<solver_type<<"  Mesh size: "<<mesh<<std::endl;
                // test(int x_mesh, int y_mesh, bool verbose);
                test unit_test(mesh, 1, false);
                // solver::solver(int order, int x_mesh, int y_mesh, int solver_type, q vector)
                solver boo (order, mesh, 1, 3, unit_test.q);
                std::vector<double> temperatures =boo.solve();
                unit_test.compare(temperatures);

                std::cout<<" 2d  " <<order<<"  order  "<<"solver_type=="<<solver_type<<"  Mesh size: "<<mesh<<std::endl;
                test new_unit_test(mesh,mesh, false);
                // solver::solver(int order, int x_mesh, int y_mesh, int solver_type, q vector)
                solver new_boo (order, mesh, mesh, solver_type, new_unit_test.q);
                std::vector<double> new_temps =new_boo.solve();
                new_unit_test.compare(new_temps);
            }
            contin:;
        }
    }
}

int main(int argc, char * argv[]){
    PetscInitialize(&argc, &argv, NULL, NULL);
    if (argc>1){
        if(std::string(argv[1])=="TEST"){make_test();}
    }
    // Declaring input variables
    int order, x_mesh, y_mesh, solver_type;
    int verification;

    PetscInitialize(&argc, &argv, NULL, NULL);
    //Calling parse to assign values to input variables
    parse(order,x_mesh,y_mesh,solver_type,verification);

    // test(int x_mesh, int y_mesh, bool verbose);
    test unit_testing(x_mesh, y_mesh, false);

    // solver::solver(int order, int x_mesh, int y_mesh, int solver_type, q vector)
    solver boo (order, x_mesh, y_mesh, solver_type, unit_testing.q);

    std::vector<double> solution= boo.solve();
    unit_testing.compare(solution);
    PetscFinalize();
    return 0;
};