#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<grvy.h>
#include<sys/time.h>
#include<time.h>
#include<unistd.h>
#include<vector>
#include<cmath>
// T(x) function
double Temperature(int xmesh, int ymesh, int i, int j);

// Q(x) function
double Q(int xmesh, int ymesh, int i, int j, double k);

// Get l2 Norm
double l2_norm_1D(const std::vector<double>& numerical, int xmesh, int ymesh);