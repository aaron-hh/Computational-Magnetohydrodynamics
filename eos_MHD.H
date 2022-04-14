#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

typedef std::array<double,10> Arrayofdouble;

class eos
{
private:
    double m_gamma;

public:
    eos();

    void set_gamma(double GAMMA);
    virtual Arrayofdouble prim_to_u(Arrayofdouble prim, double gamma, int nVar, int D);
    virtual Arrayofdouble u_to_prim(Arrayofdouble u, double gamma, int nVar, int D);
    virtual Arrayofdouble fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D);
    virtual Arrayofdouble y_fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D);
};
