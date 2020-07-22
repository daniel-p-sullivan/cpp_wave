#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "main.h"
#include "Configuration.h"
#include "DomainFunctions.h"
#include "FileWriter.h"


using std::string;      using std::cout;


int main() {
    string OUTPATH = "./out/";  //output path for solution files    
    static int count = 0;       //current solution #
    double time = 0.0;          //start time      
    double dt = 10000000;       //initialize time step for min-finding
    bool timeflag = false;      //time to stop?
    
    cell** domain;              //create the domain
    Initialize(domain);         //initialize the domain    
    InitialCondition(domain);   //set condition at t0        
    WriteFile(domain, OUTPATH, count, time);  //write the file
    
    
    for (int time_iteration = 0; time_iteration < TIME_ITER; ++time_iteration) {
        
        dt = MinRequiredTimeStep(domain, dt);    //get the time step
        if (time + dt > STOPTIME) {                 //check if time to stop
            dt = STOPTIME - time;                   //time will step onto stopping time
            timeflag = true;                        //change the flag
        }
        
        for (int rk = 0; rk < RKSTEP; ++rk) {
            UpdateGhostCells(domain, rk);               //update the GCs boundary conditions
            UpdateMainCells(domain, rk, FirstOrder);    //reconstruct the cell boundaries
            UpdateMainCells(domain, rk, ZeroFlux);      //zero out the flux
            CalculateFlux(domain, rk);                  //recalculate the flux
            UpdateMainCells(domain, rk, dt, RungeKutta);
        }
        UpdateMainCells(domain, RungeKuttaReset);   //copy final rkstep to the new bottom layer
        
        time += dt;
        if (timeflag) {
            break;
        }
        if ((time_iteration + 1) % 2 == 0) {            //lets not write 2gb of txt files
            WriteFile(domain, OUTPATH, count, time);    //write every nth time iteration
        }
    } 
    WriteFile(domain, OUTPATH, count, time);            //write the final state
    return 0;
}