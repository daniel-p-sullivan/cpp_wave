#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include "Configuration.h"
#include "main.h"

using std::string;          using std::ofstream;
using std::ostringstream;   using std::cout;
using std::endl;            using std::cerr;
using std::exception;

void WriteFile(cell** domain, string outpath, int& count, double& time) {
    
    string extension = ".txt";           //file extension to write
    
    ostringstream oss;                  //stringstream to append an int to the outpath
    oss << count << extension;          //outpath = "./out/%count.txt
    outpath += oss.str();               //setting the outpath
    const char *op = outpath.c_str();    //std::string to char*
    
    remove(op);                         //delete the previous file if it exists
    ofstream outfile(op);               //writing a new file
    
    try {                               //stay noided
   
        outfile << "XCELLS=" << XCELLS << endl          //XCELLS=%XCELLS\n
                << "YCELLS=" << YCELLS << endl          //YCELLS=%YCELLS\n
                << "GHOSTCELLS=" << GHOSTCELLS << endl  //GHOSTCELLS=%GHOSTCELLS\n
                << "TIME=" << time << endl;             //TIME=%TIME\n

        for (int i = 0; i < ROWS; ++i) {                //data input
            for (int j = 0; j < COLS; ++j) {
                cell c = domain[i][j];
                outfile << c.cx << " " << c.cy << " ";  //cx cy
                for (int k = 0; k < VARS; ++k){         //need loop through the vars at the 0th RK step
                    outfile << c.u[0][k] << " ";        //c.u[0][var1] c.u[0][var2] ... c.u[0][varn]
                }
                outfile << "\n";                        //newline
            }
        }
    } catch(exception& e) {      //thank god for c++11 b/c I have no idea what to catch here  
        cerr << e.what(); 
        cout << "fail";         //please never get here
    }
    std::cout << "Wrote " << count << extension << std::endl;
    outfile.close();    //close the file
    count++;     //increment count
};
