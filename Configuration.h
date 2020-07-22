#ifndef GUARD_CONFIGURATION_H
#define GUARD_CONFIGURATION_H
const int VARS = 3;                         //number of PDE variables
const int XCELLS = 100;                     //number of X cells
const int YCELLS = 100;                     //number of Y cells
const int GHOSTCELLS = 4;                   //number of ghost cells
static double STOPTIME = 1.0;               //stopping time [s]
static int TIME_ITER = 100000;              //number of allowed time iterations
static double XMIN = -0.1;                  //minimum x val [m]
static double YMIN = -0.1;                  //minimum y val [m]
static double XMAX = 0.1;                   //maximum x val [m]
static double YMAX = 0.1;                   //maximum y val [m]
static double COURANT = 0.1;                //Courant number
static double DX = (XMAX - XMIN) / XCELLS;  //constant DX value for structured mesh
static double DY = (YMAX - YMIN) / YCELLS;  //constant DX value for structured mesh
const int RKSTEP = 3;                       //number of Runge-Kutta steps
const int ROWS = 2 * GHOSTCELLS + YCELLS;   //calculated value for the number of rows
const int COLS = 2 * GHOSTCELLS + XCELLS;   //calculated value for the number of cols
const int NUMCELLS = ROWS * COLS;           //size of the domain, used for malloc
#endif /* CONFIGURATION_H */

