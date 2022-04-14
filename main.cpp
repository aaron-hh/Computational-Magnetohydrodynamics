//This file contains the initial settings required to run the simulations
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <array>

#include "system_MHD.H"

// Test Cases
// 0 for sod test in x direction 
// 1 for sod test in y direction
// 2 for diagonally aligned sod test
// 3 for Toro's cylindrical test
// 4 for 1-D Brio&Wu test & 2-D Brio&Wu test x-direction
// 5 2-D Brio&Wu test y-direction
// 6 for Orszag test
// 7 for Kelvin test
int test_num = 6;

// 0 for running without divergence cleaning; 1 for running with divergence cleaning
int div_clean = 1;

// Fixed at 9 for with and without divergence cleaning
double NVAR = 9;
int DIMENSION = 2;//fixed 

// simulation domain - specific to test case
int XNCELLS = 256.0;
double X0 = 0.0;
double X1 = 1.0;
int YNCELLS = 256.0;// set this to 1 when running 1-D test
double Y0 = 0.0;
double Y1 = 1.0;

// test case specific settings
double GAMMA = 5/double(3);
double TSTART = 0; 
double TEND = 0.5;

// CFL number
double C2 = 0.9;
double C1 = 0.9;

// counter for outputting DC data
int counter = 0;

int main()
{
    class system* S = new class system();
    S->setparameters(NVAR, XNCELLS, X0, X1, C1, C2, GAMMA, TSTART, TEND, DIMENSION, YNCELLS, Y0, Y1);
    S->computeInitialCondition(GAMMA, test_num);
    S->iteration(counter, test_num, div_clean);
    S->outputData(test_num);
}






