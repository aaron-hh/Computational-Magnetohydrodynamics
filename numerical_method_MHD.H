#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <array>

#include "eos_MHD.H"

typedef std::array<double,10> Arrayofdouble;

class numerical_method_MHD
{
private:
    Arrayofdouble Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);
    Arrayofdouble Minbee(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);

public:
    numerical_method_MHD();
    
    void wavespeedestimate(Arrayofdouble p, Arrayofdouble q, double n, double gamma, Arrayofdouble& wavespeed, int D);
    void ywavespeedestimate(Arrayofdouble p, Arrayofdouble q, double n, double gamma, Arrayofdouble& ywavespeed, int D);

    void compute_nplushalf_variables(Arrayofdouble a, Arrayofdouble b, Arrayofdouble c, 
    double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D);
    void ycompute_nplushalf_variables(Arrayofdouble a, Arrayofdouble b, Arrayofdouble c, 
    double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D);
    
    double compute_Bx_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch);
    double compute_By_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch);
    double compute_psix_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch);
    double compute_psiy_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch);
    Arrayofdouble compute_xflux_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch);
    Arrayofdouble compute_yflux_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch);

    void getHLLvar(Arrayofdouble a, Arrayofdouble b, Arrayofdouble p, Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& u_HLL, int nVar, int D, int div_clean);
    void ygetHLLvar(Arrayofdouble a, Arrayofdouble b, Arrayofdouble p, Arrayofdouble q, Arrayofdouble ywavespeed, Arrayofdouble& yu_HLL, int nVar, int D, int div_clean);

    void getHLLCvar_l(Arrayofdouble& u_HLL, Arrayofdouble a, Arrayofdouble p, Arrayofdouble wavespeed, Arrayofdouble& HLLCvar_l, Arrayofdouble b, Arrayofdouble q, int nVar, int D);
    void getHLLCvar_r(Arrayofdouble& u_HLL, Arrayofdouble b, Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& HLLCvar_r, Arrayofdouble a, Arrayofdouble p, int nVar, int D);
    void ygetHLLCvar_l(Arrayofdouble& yu_HLL, Arrayofdouble a, Arrayofdouble p, Arrayofdouble ywavespeed, Arrayofdouble& yHLLCvar_l, Arrayofdouble b, Arrayofdouble q, int nVar, int D);
    void ygetHLLCvar_r(Arrayofdouble& yu_HLL, Arrayofdouble b, Arrayofdouble q, Arrayofdouble ywavespeed, Arrayofdouble& yHLLCvar_r, Arrayofdouble a, Arrayofdouble p, int nVar, int D);

    void compute_HLLCflux(Arrayofdouble a, Arrayofdouble p, Arrayofdouble b, 
    Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& Fhllc, int nVar, double gamma, int D, Arrayofdouble& u_HLL, Arrayofdouble& HLLCvar_l, Arrayofdouble& HLLCvar_r, int div_clean, double Ch);
    void ycompute_HLLCflux(Arrayofdouble a, Arrayofdouble p, Arrayofdouble b, 
    Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& Fhllc, int nVar, double gamma, int D, Arrayofdouble& yu_HLL, Arrayofdouble& yHLLCvar_l, Arrayofdouble& yHLLCvar_r, int div_clean, double Ch);

    Arrayofdouble getFORCEflux(Arrayofdouble a, Arrayofdouble b,Arrayofdouble p,Arrayofdouble q, double dx, double dt, double gamma, double nVar, int D);
    Arrayofdouble ygetFORCEflux(Arrayofdouble a, Arrayofdouble b,Arrayofdouble p,Arrayofdouble q, double dx, double dt, double gamma, double nVar, int D);
};