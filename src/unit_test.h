#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#include <vector>
#include <assert.h>
#include <iostream>
#include <cmath>

class test{
    public:
        test(int x_mesh, int y_mesh, bool verbose);
        std::vector<double> q;
        std::vector<double> temps;
        void compare(std::vector<double> solution);
    private:
        double xmesh;
        double ymesh;
        double _verbose;
};
#endif