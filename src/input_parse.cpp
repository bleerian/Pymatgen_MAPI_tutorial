#include"input_parse.h"

GRVY::GRVY_Timer_Class pt; // GRVY Parsing Timer
int parse(int& order, int& x_mesh, int& y_mesh, int& solver_type, int& verification){
    pt.Init("Time Parsing");
    pt.BeginTimer("Parsing Timer");

    GRVY::GRVY_Input_Class iparse; // Input parsing object

    // Initialize/read the file
    if(! iparse.Open("./input.dat"))
    exit(1);
    // Read specific variables (with no default values provided)
    if( iparse.Read_Var("order",&order) )
    printf("--> %s = %i\n","order",order);
    if( iparse.Read_Var("x_mesh",&x_mesh) )
    printf("--> %s = %i\n","x_mesh",x_mesh);
    if( iparse.Read_Var("y_mesh",&y_mesh) )
    printf("--> %s = %i\n","y_mesh",y_mesh);
    if( iparse.Read_Var("verification",&verification) )
    printf("--> %s = %i\n","verfication mode",verification);
    if( iparse.Read_Var("solver_type",&solver_type) )
    printf("--> %s = %i\n","solver type",solver_type);
    
    // Dump the whole file to a file (appends to file if it already exists)
    printf("\n ------ Full Dump to test.out ------\n\n");
    iparse.Fdump("% ","test.out");
    printf(" ------- End Full Dump -------\n\n");
    // Close the file
    iparse.Close();
    pt.EndTimer("Parsing Timer");
    pt.Finalize();
    pt.Summarize();
    return 0;
}