#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <array>
#include <sstream>

#include "array.H"
#include "numerical_method_MHD.H"

typedef std::array<double,10> Arrayofdouble;
typedef std::vector<double> Vectorofdouble;
typedef std::vector<std::vector<double> > Vector_vectorofdouble;

class system
{
private:
    static double m_nVar;
    static int m_xnCells;
    static double m_x0;
    static double m_x1;
    static double m_dx;
    static double m_c1;
    static double m_c2;
    static double m_gamma;
    static double m_tstart;
    static double m_tend;
    static int m_d;

    static int m_ynCells;
    static double m_y0;
    static double m_y1;
    static double m_dy;

    void computeTimeStep(double m_dx, double m_dy, double m_c, double m_xnCells, double m_ynCells, double m_gamma, double& dt, double& Ch);
    void TransBC(m_array& U, double m_xnCells, double m_ynCells, int m_d);
    void PeriodicBC(m_array& U, double m_xnCells, double m_ynCells, int m_d);
    void KelvintestBC(m_array& U, double m_xnCells, double m_ynCells, int D, int m_nVar);

public:

    system(); 
    void setparameters(double NVAR, int XNCELLS, double X0, double X1, double C1, double C2, double GAMMA, double TSTART, double TEND, 
    int DIMENSION, int YNCELLS, double Y0, double Y1);
    void computeInitialCondition(double m_gamma, int test_num);
    void iteration(int& counter, int test_num, int div_clean);
    void outputData(int test_num);
};