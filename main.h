#ifndef GUARD_MAIN_H
#define GUARD_MAIN_H
#include "Configuration.h"
struct cell {    
    double dx, dy, cx, cy;        //infinitesimals dx, dy & centroid coords cx, cy
    double* left;
    double* right;
    double* top;                  //boundary values and flux
    double* bottom;
    double* flux;
    double** u;                   //PDE variables
};
int main();
#endif /* MAIN_H */

