#include "unit_test.h"

test::test(int x_mesh, int y_mesh, bool verbose): xmesh(x_mesh), ymesh(y_mesh), _verbose(verbose) {
    double k =1;
    double x;
    xmesh++;
    if (ymesh!=1){
        ymesh++;
    }
    // hard code solution and heat vector
    if(ymesh==1){
        for (int i=0; i<xmesh; i++){
            x=i;
            x=x/(xmesh-1);

            q.push_back(k*M_PI*M_PI*std::sin(M_PI*x));
            temps.push_back(k*std::sin(M_PI*x));
        }
    }else{
        for (int j=0; j<xmesh; j++){
            x=j;
            x=x/(xmesh-1);
            for (int i=0; i<ymesh; i++){
                double y = i;
                y=y/(ymesh-1);
                q.push_back(k*5*M_PI*M_PI*std::sin(M_PI*x)*std::sin(2*M_PI*y));
                temps.push_back(k*std::sin(M_PI*x)*std::sin(2*M_PI*y));
            }
        }
    }
}
void test::compare(std::vector<double> solution){
    double norm=0;
    double diff;
    assert(solution.size()==temps.size());
    for (int i = 0; i <temps.size(); ++i){
        diff= (temps[i]-solution[i])*(temps[i]-solution[i]);
        norm=norm+diff;
        if(_verbose){
            std::cout<<solution[i]<<" <----solution real----> "<<temps[i]<<std::endl;
        }
    }
    std::cout<<"L2 Error: "<<std::sqrt(norm/xmesh)<<std::endl;
}
