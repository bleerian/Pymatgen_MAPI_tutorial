#include "verification.h"

double Temperature(int xmesh, int ymesh, int x, int y) {
    if (ymesh == 1) { // 1D case
        // sol'n that fits the boundary conditions at x = 0 and x = 1
        return sin(M_PI * x / xmesh);
    } else { // 2D case
        // Asymmetrical 2D polynomial solution
        return sin(M_PI * x / xmesh) * sin(2 * M_PI * y / ymesh);
    }
}

double Q(int xmesh, int ymesh, int x, int y, double k){
    if (ymesh == 1) { // 1D case
        // sol'n that fits the boundary conditions at x = 0 and x = 1
        return k * M_PI * M_PI * sin(M_PI * x / xmesh);
    } else { // 2D case
        // Asymmetrical 2D polynomial solution
        return k * M_PI * M_PI * sin(M_PI * x / xmesh) * sin(2 * M_PI * y / ymesh);
    }
}

double l2_norm_1D(const std::vector<double>& numerical, int xmesh, int ymesh) {
    double sum = 0.0;
    int N = numerical.size();
    for (int i = 0; i < N; ++i) {
        double diff = Temperature(xmesh, ymesh, i, 0) - numerical[i];
        sum += (diff * diff);
    }
    return sqrt(sum / N);
}