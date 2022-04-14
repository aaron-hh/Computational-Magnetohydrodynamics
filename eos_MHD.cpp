//This file contains the variable conversion functions and flux functions 
#include "eos_MHD.H"

//constructor 
eos::eos()
{}

//setter for gamma 
void eos::set_gamma(double GAMMA)
{
    m_gamma = GAMMA;
}

//function to convert primitive to conservative
Arrayofdouble eos::prim_to_u(Arrayofdouble prim, double gamma, int nVar, int D)
{
  Arrayofdouble u;

  u[0] = prim[0];
  u[1] = prim[0]*prim[1];
  u[2] = prim[0]*prim[2];
  u[3] = prim[0]*prim[3];
  double v2 = prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3];
  double B2 = prim[7]*prim[7]+prim[5]*prim[5]+prim[6]*prim[6];
  u[4] = (prim[4]/(gamma-1))+(0.5*prim[0]*v2)+(0.5*B2);
  u[5] = prim[5];
  u[6] = prim[6];
  u[7] = prim[7];

  u[8] = prim[8];

  return u;
} 

//function to convert conservative to primitive
Arrayofdouble eos::u_to_prim(Arrayofdouble u, double gamma, int nVar, int D)
{  
  Arrayofdouble prim;

  prim[0] = u[0];//rho
  prim[1] = u[1]/u[0];//vx
  prim[2] = u[2]/u[0];//vy
  prim[3] = u[3]/u[0];//vz
  double v2 = prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3];
  double B2 = u[7]*u[7]+u[5]*u[5]+u[6]*u[6];
  prim[4] = (u[4]-(0.5*u[0]*v2)-(0.5*B2))*(gamma-1);//p
  prim[5] = u[5];//By
  prim[6] = u[6];//Bz
  prim[7] = u[7];//Bx

  prim[8] = u[8];

  return prim;
}

//defining fluxf function
Arrayofdouble eos::fluxf(Arrayofdouble u,Arrayofdouble prim, int nVar, int D)
{
  Arrayofdouble fluxfunction;

  double B2 = prim[7]*prim[7]+prim[5]*prim[5]+prim[6]*prim[6];

  fluxfunction[0] = prim[0]*prim[1];
  fluxfunction[1] = (prim[0]*prim[1]*prim[1])+prim[4]+(0.5*B2)-(prim[7]*prim[7]);
  fluxfunction[2] = prim[0]*prim[1]*prim[2] - prim[7]*prim[5];
  fluxfunction[3] = prim[0]*prim[1]*prim[3] - prim[7]*prim[6];
  fluxfunction[4] = (u[4]+prim[4]+0.5*B2)*prim[1] - (prim[1]*prim[7] + prim[2]*prim[5] + prim[3]*prim[6])*prim[7];
  fluxfunction[5] = prim[5]*prim[1] - prim[7]*prim[2];
  fluxfunction[6] = prim[6]*prim[1] - prim[7]*prim[3];
  fluxfunction[7] = 0;

  fluxfunction[8] = 0;

  return fluxfunction;
}

//define yflux function
Arrayofdouble eos::y_fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D)
{
  Arrayofdouble yfluxfunction;

  double B2 = prim[7]*prim[7]+prim[5]*prim[5]+prim[6]*prim[6];

  yfluxfunction[0] = prim[0]*prim[2];
  yfluxfunction[1] = prim[0]*prim[1]*prim[2] - prim[7]*prim[5];
  yfluxfunction[2] = (prim[0]*prim[2]*prim[2]) + prim[4] + (0.5*B2) - (prim[5]*prim[5]);
  yfluxfunction[3] = prim[0]*prim[2]*prim[3] - prim[5]*prim[6];
  yfluxfunction[4] = (u[4]+prim[4]+0.5*B2)*prim[2] - (prim[1]*prim[7] + prim[2]*prim[5] + prim[3]*prim[6])*prim[5];
  yfluxfunction[5] = 0;
  yfluxfunction[6] = prim[6]*prim[2] - prim[5]*prim[3];
  yfluxfunction[7] = prim[7]*prim[2] - prim[5]*prim[1];

  yfluxfunction[8] = 0;

  return yfluxfunction;
}
