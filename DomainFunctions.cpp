#include "main.h"
#include "DomainFunctions.h"
#include "Configuration.h"
#include "math.h"
#include <iostream>
#include <tgmath.h>
#include <limits>
#include <cmath>
#include <stdlib.h>
using std::abs;     using std::max;
using std::cout;    
void UpdateDomain(cell**& domain, void function(cell&)) {         //call function on entire domain
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            function(domain[i][j]);
        }
    }
}
void UpdateGhostCells(cell**& domain, int& rk) {     //call function just on the ghost cells
    
    for (int i = 0; i < COLS; ++i) {                //left and right
        for (int gc = 0; gc < GHOSTCELLS; ++gc) {
            //function(domain[i][j], rk);
            int variable = 0;                       //set up w/o a loop so you can change the
            domain[gc][i].u[rk][variable]           //boundary condition based on the variable
                    = domain[2 * GHOSTCELLS - gc - 1][i].u[rk][variable];
            
            domain[YCELLS + GHOSTCELLS + gc][i].u[rk][variable]
                    = domain[YCELLS + GHOSTCELLS - gc - 1][i].u[rk][variable];
            
            variable = 1;
            domain[gc][i].u[rk][variable]
                    = -domain[2 * GHOSTCELLS - gc - 1][i].u[rk][variable];
            
            domain[YCELLS + GHOSTCELLS + gc][i].u[rk][variable]
                    = -domain[YCELLS + GHOSTCELLS - gc - 1][i].u[rk][variable];
            
            variable = 2;
            domain[gc][i].u[rk][variable]
                    = domain[2 * GHOSTCELLS - gc - 1][i].u[rk][variable];
            
            domain[YCELLS + GHOSTCELLS + gc][i].u[rk][variable]
                    = domain[YCELLS + GHOSTCELLS - gc - 1][i].u[rk][variable];           
        }
    }
    //top and bottom
    for (int i = 0; i < ROWS; ++i) {
        for (int gc = 0; gc < GHOSTCELLS; ++gc) {
            //function(domain[i][j], rk);
            int variable = 0;                       //set up w/o a loop so you can change the
            domain[i][gc].u[rk][variable]           //boundary condition based on the variable
                    = domain[i][2 * GHOSTCELLS - gc - 1].u[rk][variable];
            
            domain[i][XCELLS + GHOSTCELLS + gc].u[rk][variable]
                    = domain[i][XCELLS + GHOSTCELLS - gc - 1].u[rk][variable];
            
            variable = 1;
            domain[i][gc].u[rk][variable]
                    = -domain[i][2 * GHOSTCELLS - gc - 1].u[rk][variable];
            
            domain[i][XCELLS + GHOSTCELLS + gc].u[rk][variable]
                    = -domain[i][XCELLS + GHOSTCELLS - gc - 1].u[rk][variable];
            
            variable = 2;
            domain[i][gc].u[rk][variable]
                    = -domain[i][2 * GHOSTCELLS - gc - 1].u[rk][variable];
            
            domain[i][XCELLS + GHOSTCELLS + gc].u[rk][variable]
                    = -domain[i][XCELLS + GHOSTCELLS - gc - 1].u[rk][variable];
        }
    }
    
}
void UpdateMainCells(cell**& domain, int& rk, void function(cell&, int&)) {      //call function on the main cells
    for (int i = GHOSTCELLS - 1; i < ROWS - GHOSTCELLS + 1; ++i) {
        for (int j = GHOSTCELLS - 1; j < COLS - GHOSTCELLS + 1; ++j) {
            function(domain[i][j], rk);
        }
    }
}
void UpdateMainCells(cell**& domain, int& rk, double& dt, void function(cell&, int&, double&)) {      //call function on the main cells
    for (int i = GHOSTCELLS; i < ROWS - GHOSTCELLS; ++i) {
        for (int j = GHOSTCELLS; j < COLS - GHOSTCELLS; ++j) {
            function(domain[i][j], rk, dt);
        }
    }
}
void UpdateMainCells(cell**& domain, void function(cell&)) {
    for (int i = GHOSTCELLS; i < ROWS - GHOSTCELLS; ++i) {
        for (int j = GHOSTCELLS; j < COLS - GHOSTCELLS; ++j) {
            function(domain[i][j]);
        }
    }
}
void CalculateFlux(cell** domain, int& rk) {
    
    double dx = domain[0][0].dx;    //structured => dx, dy are const
    double dy = domain[0][0].dy;    //leaving here for implementation of unstructured
    
    for (int j = GHOSTCELLS; j < XCELLS + GHOSTCELLS; ++j) {            //vertical interface discrete flux calc
        for (int interface = GHOSTCELLS; 
                 interface < YCELLS + GHOSTCELLS + 1; ++interface) {
            double* lv = domain[interface - 1][j].right;
            double* rv = domain[interface][j].left;
            double g = 9.81;
            double hl = lv[0];
            double ul = lv[1] / hl;
            double cl = sqrt(g * hl);
            double vmaxl = abs(ul) + cl;
            
            double hr = rv[0];
            double ur = rv[1] / hr;
            double cr = sqrt(g * hr);
            double vmaxr = abs(ur) + cr;
            double vmax = max(vmaxr, vmaxl);
            
            double* FL = F(lv);
            double* FR = F(rv);
            
            for(int variable = 0; variable < VARS; variable++) {
                double f = 0.5 * (FL[variable] + FR[variable])
                           - 0.5 * vmax * (rv[variable] - lv[variable]);
                f *= dy;
                domain[interface - 1][j].flux[variable] -= f;
                domain[interface][j].flux[variable] += f;
            }
            delete[] FL;
            delete[] FR;
        }
    }
    
    for (int i = GHOSTCELLS; i < YCELLS + GHOSTCELLS; ++i) {            //horizontal interface discrete flux calc
        for (int interface = GHOSTCELLS; 
                 interface < XCELLS + GHOSTCELLS + 1; ++interface) {
            double* bv = domain[i][interface - 1].top;
            double* tv = domain[i][interface].bottom;
            double g = 9.81;
            double hb = bv[0];
            double vb = bv[2] / hb;
            double cb = sqrt(g * hb);
            double vmaxb = abs(vb) + cb;
            
            double ht = tv[0];
            double ut = tv[1] / ht;
            double ct = sqrt(g * ht);
            double vmaxt = abs(ut) + ct;
            double vmax = max(vmaxt, vmaxb);
            
            double* GB = G(bv);
            double* GT = G(tv);
            
            for(int variable = 0; variable < VARS; variable++) {
                double f = 0.5 * (GB[variable] + GT[variable])
                           - 0.5 * vmax * (tv[variable] - bv[variable]);
                f *= dx;
                domain[i][interface - 1].flux[variable] -= f;
                domain[i][interface].flux[variable] += f;
            }
            delete[] GB;
            delete[] GT;
        }
    }    
}
void FirstOrder(cell& c, int& rk) {                        //First Order reconstruction method
    for (int variable = 0; variable < VARS; ++variable) {   //doesn't attempt to predict the values
        c.left[variable] = c.u[rk][variable];               //at the boundary of the cell. /* SEGMENTATION FAULT HERE */
        c.right[variable] = c.u[rk][variable];              //sets the boundary value to
        c.top[variable] = c.u[rk][variable];                //the value at the centroid.
        c.bottom[variable] = c.u[rk][variable];
    }
}
void WENO(cell& c, int& rk) {
    
    
}
void RungeKutta(cell& c, int& rk, double& dt) {
    for (int variable = 0; variable < VARS; ++variable) {
        double da = c.dx * c.dy;
        switch(rk) {
            case 0:
                c.u[rk + 1][variable] = c.u[rk][variable] + dt / da * c.flux[variable];
                break;
            case 1:
                c.u[rk + 1][variable] = 0.75 * c.u[rk - 1][variable] + 
                                        0.25 * (c.u[rk][variable] + dt / da * c.flux[variable]);
                break;
            case 2:
                c.u[rk + 1][variable] = 1.0 / 3.0 * c.u[rk - 2][variable] + 
                                        2.0 / 3.0 * (c.u[rk][variable] + dt / da * c.flux[variable]);
                break;
            default:
                cout << "broken RK";
        }    
    }   
}
void RungeKuttaReset(cell& c) {
    for (int variable = 0; variable < VARS; ++variable) {
        c.u[0][variable] = c.u[RKSTEP][variable];
    }  
}
void ZeroFlux(cell& c, int& rk) {
    for (int variable = 0; variable < VARS; ++variable) {
        c.flux[variable] = 0.0;
    }
}
void Initialize(cell**& domain) {
        
    domain = new cell*[ROWS];
    for (int h = 0; h < ROWS; ++h) {
        domain[h] = new cell[COLS];
    }

    
    //set initial dx, dy, cx, cy values for all cells. ONLY WORKS WITH STRUCTURED DOMAIN
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            domain[i][j].dx = DX;
            domain[i][j].dy = DY;
            domain[i][j].cx = XMIN + DX * (j + 0.5 - GHOSTCELLS);  //indexing formula for the centroids
            domain[i][j].cy = YMIN + DY * (i + 0.5 - GHOSTCELLS);
            domain[i][j].left = new double[VARS];                  //allocate mem for the cell double*
            domain[i][j].right = new double[VARS];
            domain[i][j].top = new double[VARS];
            domain[i][j].bottom = new double[VARS];
            domain[i][j].flux = new double[VARS];
            domain[i][j].u = new double*[RKSTEP+1];                  //make the double**
            for (int k = 0; k <= RKSTEP; ++k) {                     //easier to access later with u[][] syntax
                domain[i][j].u[k] = new double[VARS];              //have to allocate a double** w/i a cell**
                for (int l = 0; l < VARS; ++l) {
                    domain[i][j].u[k][l] = 1E-15;         //assign values w/i the double**
                }
            }
            for (int m = 0; m < VARS; ++m){                    //assign values to all the double*
                domain[i][j].left[m] = 1E-15;             //will be overwritten
                domain[i][j].right[m] = 1E-15;            //in InitialCondition
                domain[i][j].top[m] = 1E-15;
                domain[i][j].bottom[m] = 1E-15;
                domain[i][j].flux[m] = 1E-15;              
            }
        }
    }    
}
void InitialCondition(cell**& domain) {
    
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j){
            for (int rk = 0; rk <= RKSTEP; ++rk) {
                double r = domain[i][j].cx * domain[i][j].cx + domain[i][j].cy * domain[i][j].cy;                 //radius from the center
                double amplitude = 3.0;
                //double letsseewhatthisdoes = 0.001;
                double variance = 0.2 * 0.2;
                
                domain[i][j].u[rk][0] = amplitude * exp(-r / variance);
                domain[i][j].u[rk][1] = 0.0;//letsseewhatthisdoes * -sqrt(r) * sin(atan2(domain[i][j].cy, domain[i][j].cx));
                domain[i][j].u[rk][2] = 0.0;//letsseewhatthisdoes * sqrt(r) * cos(atan2(domain[i][j].cy, domain[i][j].cx));  
            }                           
        }
    }    
}
void Delete(cell& c) {
    delete[] c.left;
    delete[] c.right;
    delete[] c.top;
    delete[] c.bottom;
    delete[] c.flux;
    for (int k = 0; k <= RKSTEP; ++k) {
        delete[] c.u[k];
    }
    delete[] c.u;
}
double& MinRequiredTimeStep(cell**& domain, double& dt) { //what is the minimum required time step
                                             
    dt = CellTimeStep(domain[0][0]);
    for (int i = GHOSTCELLS; i < ROWS - GHOSTCELLS; ++i) {
        for (int j = GHOSTCELLS; j < COLS - GHOSTCELLS; ++j) {
            double c_dt = CellTimeStep(domain[i][j]);
            if (c_dt < dt) {
                dt = c_dt;
            }    
        }
    }
    return dt;
}
double CellTimeStep(cell& c) {
    return 1.0 / (1.0 / (c.dx / (abs(c.u[0][1] / c.u[0][0]) + sqrt(9.81 * c.u[0][0])))
                + 1.0 / (c.dy / (abs(c.u[0][2] / c.u[0][0]) + sqrt(9.81 * c.u[0][0])))); //im a wizard
}
double* F(double*& x) { //From Navier-Stokes, conservation of mass, incompressible
    double* ret = new double[3];
    double h = x[0];
    double u = x[1] / h;
    double v = x[2] / h;
    ret[0] = h * u;
    ret[1] = h * u * u + 0.5 * 9.81 * h * h;
    ret[2] = h * u * v;
    return ret;
}
double* G(double*& x) { //From Navier-Stokes, conservation of mass, incompressible
    double* ret = new double[3];
    double h = x[0];
    double u = x[1] / h;
    double v = x[2] / h;
    ret[0] = h * v;
    ret[1] = h * u * v;
    ret[2] = h * v * v + 0.5 * 9.81 * h * h;
    return ret;
}