//This file contains the iteration loop, compute time step function, boundary conditions and the initial data for all the MHD test cases
#include "system_MHD.H"

using namespace std;

//initiating classes
numerical_method_MHD* NM = new numerical_method_MHD();
eos* E = new eos();

//data storage for 1D
m_array U;
m_array PRIM;
m_array UPLUS1;
m_array UIL_NPLUSHALF;
m_array UIR_NPLUSHALF;
m_array PRIML_NPLUSHALF;
m_array PRIMR_NPLUSHALF;
m_array WAVESPEED;
m_array FLUX;

//extra data storage for 2D
m_array UBAR;
m_array PRIMBAR;
m_array UIL_NPLUSHALF_BAR;
m_array UIR_NPLUSHALF_BAR;
m_array PRIML_NPLUSHALF_BAR;
m_array PRIMR_NPLUSHALF_BAR;
m_array WAVESPEED_BAR;
m_array YFLUX;

//extra data storage for MHD
m_array U_HLL;
m_array HLLCVAR_L;
m_array HLLCVAR_R;
m_array YU_HLL;
m_array YHLLCVAR_L;
m_array YHLLCVAR_R;
m_array UBAR2;

//data storage for computing time step
m_array CS_2;
m_array CA_2X;
m_array CF_X;
m_array CF_Y;
m_array CF_Z;

//allocating memory space for static
double system::m_nVar;
int system::m_xnCells;
double system::m_x0;
double system::m_x1;
double system::m_dx;
double system::m_c1;
double system::m_c2;
double system::m_gamma;
double system::m_tstart;
double system::m_tend;
int system::m_d;

//2D
int system::m_ynCells;
double system::m_y0;
double system::m_y1;
double system::m_dy;

//constructor
system::system()
{}

//define a function for parameter settings
void system::setparameters(double NVAR, int XNCELLS, double X0, double X1,double C1, double C2, double GAMMA, double TSTART, 
double TEND, int DIMENSION, int YNCELLS, double Y0, double Y1)
{
    m_nVar = NVAR;
    m_xnCells = XNCELLS;
    m_x0 = X0;
    m_x1 = X1;
    m_dx = (X1-X0)/XNCELLS;
    m_c1 = C1;
    m_c2 = C2;
    m_gamma = GAMMA;
    m_tstart = TSTART;
    m_tend = TEND; 
    m_d = DIMENSION;

    //2D
    m_ynCells = YNCELLS;
    m_y0 = Y0;
    m_y1 = Y1;
    m_dy = (Y1-Y0)/YNCELLS;

    U.setSize(m_xnCells+4, m_ynCells+4, 1);
    PRIM.setSize(m_xnCells+4,  m_ynCells+4, 1);
    UPLUS1.setSize(m_xnCells+4,  m_ynCells+4, 1);
    UIL_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1);
    UIR_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1);
    PRIML_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1);
    PRIMR_NPLUSHALF.setSize(m_xnCells+2,  m_ynCells+2, 1);
    WAVESPEED.setSize(m_xnCells+1,  m_ynCells+1, 1);
    FLUX.setSize(m_xnCells+1,  m_ynCells+1, 1);

    //2D
    UBAR.setSize(m_xnCells+4,  m_ynCells+4, 1);
    PRIMBAR.setSize(m_xnCells+4,  m_ynCells+4, 1);
    UIL_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1);
    UIR_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1);
    PRIML_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1);
    PRIMR_NPLUSHALF_BAR.setSize(m_xnCells+2,  m_ynCells+2, 1);
    WAVESPEED_BAR.setSize(m_xnCells+1,  m_ynCells+1, 1);
    YFLUX.setSize(m_xnCells+1,  m_ynCells+1, 1);

    //MHD
    U_HLL.setSize(m_xnCells+2,  m_ynCells+2, 1);
    HLLCVAR_L.setSize(m_xnCells+2,  m_ynCells+2, 1);
    HLLCVAR_R.setSize(m_xnCells+2,  m_ynCells+2, 1);
    YU_HLL.setSize(m_xnCells+2,  m_ynCells+2, 1);
    YHLLCVAR_L.setSize(m_xnCells+2,  m_ynCells+2, 1);
    YHLLCVAR_R.setSize(m_xnCells+2,  m_ynCells+2, 1);
    UBAR2.setSize(m_xnCells+4,  m_ynCells+4, 1);

    //time step
    CS_2.setSize(m_xnCells+4, m_ynCells+4, 1);
    CA_2X.setSize(m_xnCells+4, m_ynCells+4, 1);
    CF_X.setSize(m_xnCells+4, m_ynCells+4, 1);
    CF_Y.setSize(m_xnCells+4, m_ynCells+4, 1);
    CF_Z.setSize(m_xnCells+4, m_ynCells+4, 1);
}

//compute time step
void system::computeTimeStep(double m_dx, double m_dy, double m_c, double m_xnCells, double m_ynCells, double m_gamma, double& dt, double& Ch)
{
  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
    //computing Ca
    CA_2X(i,j,0,0) = (PRIM(i,j,0,7)*PRIM(i,j,0,7) + PRIM(i,j,0,5)*PRIM(i,j,0,5) + PRIM(i,j,0,6)*PRIM(i,j,0,6)) / (PRIM(i,j,0,0));

    //computing Cs
    CS_2(i,j,0,0) = PRIM(i,j,0,4)*m_gamma/PRIM(i,j,0,0);

    //computing Cf
    CF_X(i,j,0,0) = sqrt(0.5*((CS_2(i,j,0,0)+CA_2X(i,j,0,0)) + sqrt(pow((CS_2(i,j,0,0)+CA_2X(i,j,0,0)),2) - (4*CS_2(i,j,0,0)*pow(PRIM(i,j,0,7),2)/PRIM(i,j,0,0)))));
    CF_Y(i,j,0,0) = sqrt(0.5*((CS_2(i,j,0,0)+CA_2X(i,j,0,0)) + sqrt(pow((CS_2(i,j,0,0)+CA_2X(i,j,0,0)),2) - (4*CS_2(i,j,0,0)*pow(PRIM(i,j,0,5),2)/PRIM(i,j,0,0)))));
    }
  }

  //finding Ch by getting the maximum between Ch in x, y ,z directions
  double Ch1 = 0;
  double Ch_x = 0;
  double Ch_y = 0;
  double Ch_z = 0;

  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
      if((CF_X(i,j,0,0)+fabs(PRIM(i,j,0,1)))>Ch_x)
      {
      Ch_x = CF_X(i,j,0,0)+fabs(PRIM(i,j,0,1));
      }
      if((CF_Y(i,j,0,0)+fabs(PRIM(i,j,0,2)))>Ch_y)
      {
      Ch_y = CF_Y(i,j,0,0)+fabs(PRIM(i,j,0,2));
      }
    }
  }
    
  Ch1 = std::max(Ch_x, Ch_y);
  Ch = std::max(Ch1, Ch_z);

  //computing time step for iteration
  dt = m_c*std::min(m_dx,m_dy)/Ch;
}

//defining initial conditions in primitive variable form
void system::computeInitialCondition(double m_gamma, int test_num)
{
  Arrayofdouble IC_L, IC_R;
  double x, y;

  switch(test_num)
  {
    case(0):
    {
    //Sod test x direction
    Arrayofdouble IC_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    Arrayofdouble IC_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};

    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        if(x<0.5)
        {
          PRIM(i,j,0) = IC_L;
        }
        else
        {
          PRIM(i,j,0) = IC_R;
        }
      }
    }
    break;
    }

    case(1):
    {
    //Sod test y direction
    Arrayofdouble IC_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    Arrayofdouble IC_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};

    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        if(y<0.5)
        {
          PRIM(i,j,0) = IC_L;
        }
        else
        {
          PRIM(i,j,0) = IC_R;
        }
      }
    }
    break;
    }

    case(2):
    {
    //Diagonally aligned sod test
    Arrayofdouble IC_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    Arrayofdouble IC_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};

    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        if(y<(1-x))
        {
          PRIM(i,j,0) = IC_L;
        }
        else
        {
          PRIM(i,j,0) = IC_R;
        }
      }
    }
    break;
    }

    case(3):
    {
    //Toro's cylindrical test case
    Arrayofdouble IC_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    Arrayofdouble IC_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};

    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        if(sqrt(pow((x-1),2)+pow((y-1),2))<0.4)
        {
          PRIM(i,j,0) = IC_L;
        }
        else
        {
          PRIM(i,j,0) = IC_R;
        }
      }
    }
    break;
    }

    case(4):
    {
    //Brio & Wu test x-direction
    Arrayofdouble IC_L = {1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.75};
    Arrayofdouble IC_R = {0.125, 0.0, 0.0, 0.0, 0.1, -1.0, 0.0, 0.75};

    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        if(x<400)
        {
          PRIM(i,j,0) = IC_L;
        }
        else
        {
          PRIM(i,j,0) = IC_R;
        }
      }
    }

    break;
    }

    case(5):
    {
    //Brio & Wu test y-direction
    Arrayofdouble IC_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 0.0, 1.0};
    Arrayofdouble IC_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.75, 0.0, -1.0};

    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        if(y<400)
        {
          PRIM(i,j,0) = IC_L;
        }
        else
        {
          PRIM(i,j,0) = IC_R;
        }
      }
    }

    break;
    }

    case(6):
    {
    //Orszag-tang test
    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        Arrayofdouble IC = {m_gamma*m_gamma, -sin(2*M_PI*y), sin(2*M_PI*x), 0.0, m_gamma, sin(4*M_PI*x), 0.0, -sin(2*M_PI*y),0.0};
        PRIM(i,j,0) = IC;
      }
    }
    break;
    }

    case(7):
    {
    //Kelvin test
    for(int i=0; i<m_xnCells+4; i++)
    {
      for(int j=0; j<m_ynCells+4; j++)
      {
        double x = m_x0 + (i-1)*m_dx;
        double y = m_y0 + (j-1)*m_dy;

        Arrayofdouble IC = {1.0, 0.5*(tanh(20*y)), 0.1*sin(2*M_PI*x)*(exp(-pow(y,2)/0.01)), 0.0, 1/m_gamma, 0.0, 0.1*sin(M_PI/3), 0.1*cos(M_PI/3), 0.0};
        PRIM(i,j,0) = IC;
      }
    }
    break;
    }
  }

}

//defining function for applying TRANSMISSIVE boundary condition 
void system::TransBC(m_array& U, double m_xnCells, double m_ynCells, int D)
{   
  //transmissive left right
  for(int j=0; j<m_ynCells+4; j++) 
  {
  U(0,j,0) = U(2,j,0);
  U(1,j,0) = U(2,j,0);
  U(m_xnCells+2,j,0) = U(m_xnCells+1,j,0);
  U(m_xnCells+3,j,0) = U(m_xnCells+1,j,0);
  }

  if(D == 2)
  {
  //transmissive top bottom
  for(int i=0; i<m_xnCells+4; i++)
  {
  U(i,0,0) = U(i,2,0);
  U(i,1,0) = U(i,2,0);
  U(i,m_ynCells+2,0) = U(i,m_ynCells+1,0);
  U(i,m_ynCells+3,0) = U(i,m_ynCells+1,0);
  }
  }
}

//defining function for applying PERIODIC boundary condition 
void system::PeriodicBC(m_array& U, double m_xnCells, double m_ynCells, int D)
{   
//periodic boundary condition L-R
  for(int j=0; j<m_ynCells+4; j++) 
  {
  U(0,j,0) = U(m_xnCells,j,0);
  U(1,j,0) = U(m_xnCells+1,j,0);
  U(m_xnCells+3,j,0) = U(3,j,0);
  U(m_xnCells+2,j,0) = U(2,j,0);
  }

  //periodic boundary condition T-B
  for(int i=0; i<m_xnCells+4; i++) 
  {
  U(i,0,0) = U(i,m_xnCells,0);
  U(i,1,0) = U(i,m_xnCells+1,0);
  U(i,m_xnCells+3,0) = U(i,3,0);
  U(i,m_xnCells+2,0) = U(i,2,0);
  }
}

//defining function for applying PERIODIC boundary condition 
void system::KelvintestBC(m_array& U, double m_xnCells, double m_ynCells, int D, int m_nVar)
{   
//periodic boundary condition L-R
  for(int j=0; j<m_ynCells+4; j++) 
  {
  U(0,j,0) = U(m_xnCells,j,0);
  U(1,j,0) = U(m_xnCells+1,j,0);
  U(m_xnCells+3,j,0) = U(3,j,0);
  U(m_xnCells+2,j,0) = U(2,j,0);
  }

  //reflective boundary condition T-B
  for(int i=0; i<m_xnCells+4; i++) 
  {
  //ghost vy = - vy
  //ghost By = - By
  //velocity y
  U(i,0,0,2) = -U(i,3,0,2);
  U(i,1,0,2) = -U(i,2,0,2);
  U(i,m_ynCells+2,0,2) = -U(i,m_ynCells+1,0,2);
  U(i,m_ynCells+3,0,2) = -U(i,m_ynCells,0,2);

  //By
  U(i,0,0,5) = -U(i,3,0,5);
  U(i,1,0,5) = -U(i,2,0,5);
  U(i,m_ynCells+2,0,5) = -U(i,m_ynCells+1,0,5);
  U(i,m_ynCells+3,0,5) = -U(i,m_ynCells,0,5);

  //rho, vx, vz, p, bx, bz
  for(int l=0; l<m_nVar; l++)
  {
  if(m_nVar!=2 && m_nVar!=5)
  {
  U(i,0,0,l) = U(i,3,0,l);
  U(i,1,0,l) = U(i,2,0,l);
  U(i,m_ynCells+2,0,l) = U(i,m_ynCells+1,0,l);
  U(i,m_ynCells+3,0,l) = U(i,m_ynCells,0,l);
  }
  }
  }
}

//defining function for iteration
void system::iteration(int& counter, int test_num, int div_clean)
{
Vector_vectorofdouble dc;
double t = m_tstart;
double dx = m_dx;
double dy = m_dy;
double Ch = 0.0;
double dt;
E->set_gamma(m_gamma);

//convert from primitive to conservative
for(int i=0; i<m_xnCells+4; i++)
{
  for(int j=0; j<m_ynCells+4; j++)
  {
    U(i,j,0) = E->prim_to_u(PRIM(i,j,0), m_gamma, m_nVar, m_d);
  }
}

do
{
  //use lower CFL number for the first few iteration for stability
  counter +=1;
  double C;
  if(counter<20)
  {
    C = m_c1;
  }
  else
  {
    C = m_c2;
  }

  //compute the stable time step for this iteration
  computeTimeStep(m_dx, m_dy, C, m_xnCells, m_ynCells, m_gamma, dt, Ch);
  t = t + dt;
  std::cout<<t<<std::endl;

  //applying boundary condition
  if(test_num == 6)//Orszag test case
  {
    PeriodicBC(U, m_xnCells, m_ynCells, m_d);
  }
  else if(test_num == 7)//Kelvin test case
  {
    KelvintestBC(U, m_xnCells, m_ynCells, m_d, m_nVar);
  }
  else
  {
    TransBC(U, m_xnCells, m_ynCells, m_d);
  }

  //calculate uiL_nplushalf, uiR_nplushalf for x cells (MUSCL-Hancock)
   for(int i = 0; i<m_xnCells+2; i++)
   {
     for(int j = 0; j<m_ynCells+2; j++)
     {
       NM->compute_nplushalf_variables(U(i,j+2,0), U(i+1,j+2,0), U(i+2,j+2,0), m_dx, dt, m_gamma, UIL_NPLUSHALF(i,j,0), UIR_NPLUSHALF(i,j,0), m_nVar, m_d);

       PRIML_NPLUSHALF(i,j,0) = E->u_to_prim(UIL_NPLUSHALF(i,j,0),m_gamma,m_nVar,m_d);
       PRIMR_NPLUSHALF(i,j,0) = E->u_to_prim(UIR_NPLUSHALF(i,j,0),m_gamma,m_nVar,m_d);
     }
  }

  //calculate xHLLCflux
  for(int i = 0; i<m_xnCells+1; i++)
  {
    for(int j = 0; j<m_ynCells+1; j++)
    {
    //HLLC solver
    NM->wavespeedestimate(PRIMR_NPLUSHALF(i,j,0), PRIML_NPLUSHALF(i+1,j,0), m_xnCells, m_gamma, WAVESPEED(i,j,0), m_d);
    NM->compute_HLLCflux(UIR_NPLUSHALF(i,j,0), PRIMR_NPLUSHALF(i,j,0), UIL_NPLUSHALF(i+1,j,0), 
    PRIML_NPLUSHALF(i+1,j,0), WAVESPEED(i,j,0), FLUX(i,j,0), m_nVar, m_gamma, m_d, U_HLL(i,j,0), HLLCVAR_L(i,j,0), HLLCVAR_R(i,j,0), div_clean, Ch);
    
    //FORCE solver
    //FLUX(i,j,0) = NM->getFORCEflux(UIR_NPLUSHALF(i,j,0), UIL_NPLUSHALF(i+1,j,0), PRIMR_NPLUSHALF(i,j,0), PRIML_NPLUSHALF(i+1,j,0), m_dx, dt, m_gamma, m_nVar, m_d);
    }
  }

  // x update the data using Godunov
  for(int i = 2; i<m_xnCells+2; i++)
  {
    for(int j = 2; j<m_ynCells+2; j++)
    {
      for(int k = 0; k<m_nVar; k++)
      {
      UBAR(i,j,0,k) = U(i,j,0,k) - (dt/m_dx) * (FLUX(i-1,j-2,0,k) - FLUX(i-2,j-2,0,k));
      }
    }
  }

  //applying boundary condition
  if(test_num == 6)//Orszag test case
  {
    PeriodicBC(UBAR, m_xnCells, m_ynCells, m_d);
  }
  else if(test_num == 7)//Kelvin test case
  {
    KelvintestBC(UBAR, m_xnCells, m_ynCells, m_d, m_nVar);
  }
  else
  {
    TransBC(UBAR, m_xnCells, m_ynCells, m_d);
  }

  //calculate uiL_nplushalf, uiR_nplushalf for y cells (MUSCL-Hancock)
   for(int i = 0; i<m_xnCells+2; i++)
   {
     for(int j = 0; j<m_ynCells+2; j++)
     {
       NM->ycompute_nplushalf_variables(UBAR(i+2,j,0), UBAR(i+2,j+1,0), UBAR(i+2,j+2,0), m_dy, dt, m_gamma, UIL_NPLUSHALF_BAR(i,j,0), UIR_NPLUSHALF_BAR(i,j,0), m_nVar, m_d);

       PRIML_NPLUSHALF_BAR(i,j,0) = E->u_to_prim(UIL_NPLUSHALF_BAR(i,j,0),m_gamma,m_nVar,m_d);
       PRIMR_NPLUSHALF_BAR(i,j,0) = E->u_to_prim(UIR_NPLUSHALF_BAR(i,j,0),m_gamma,m_nVar,m_d);
     }
  }

  //calculate yHLLCflux
  for(int i = 0; i<m_xnCells+1; i++)
  {
    for(int j = 0; j<m_ynCells+1; j++)
    {
    //HLLC solver
    NM->ywavespeedestimate(PRIMR_NPLUSHALF_BAR(i,j,0), PRIML_NPLUSHALF_BAR(i,j+1,0), m_ynCells, m_gamma, WAVESPEED_BAR(i,j,0), m_d);
    NM->ycompute_HLLCflux(UIR_NPLUSHALF_BAR(i,j,0), PRIMR_NPLUSHALF_BAR(i,j,0), UIL_NPLUSHALF_BAR(i,j+1,0), 
    PRIML_NPLUSHALF_BAR(i,j+1,0), WAVESPEED_BAR(i,j,0), YFLUX(i,j,0), m_nVar, m_gamma, m_d, YU_HLL(i,j,0), YHLLCVAR_L(i,j,0), YHLLCVAR_R(i,j,0), div_clean, Ch);
    
    //FORCE solver
    //YFLUX(i,j,0) = NM->ygetFORCEflux(UIR_NPLUSHALF_BAR(i,j,0), UIL_NPLUSHALF_BAR(i,j+1,0), PRIMR_NPLUSHALF_BAR(i,j,0), PRIML_NPLUSHALF_BAR(i,j+1,0), m_dy, dt, m_gamma, m_nVar, m_d);
    
    }
  }

  if(div_clean == 0)// no divergence cleaning
  {
    // y update the data using Godunov
    for(int i = 2; i<m_xnCells+2; i++)
    {
      for(int j = 2; j<m_ynCells+2; j++)
      {
        for(int k = 0; k<m_nVar; k++)
        {
        U(i,j,0,k) = UBAR(i,j,0,k) - (dt/m_dy) * (YFLUX(i-2,j-1,0,k) - YFLUX(i-2,j-2,0,k));
        }
      }
    }
  }

  else if(div_clean == 1)// divergence cleaning
  {
    // y update the data using Godunov
    for(int i = 2; i<m_xnCells+2; i++)
    {
      for(int j = 2; j<m_ynCells+2; j++)
      {
        for(int k = 0; k<m_nVar; k++)
        {
        UBAR2(i,j,0,k) = UBAR(i,j,0,k) - (dt/m_dy) * (YFLUX(i-2,j-1,0,k) - YFLUX(i-2,j-2,0,k));
        }
      }
    }

    //applying boundary condition
    if(test_num == 6)//Orszag test case
    {
      PeriodicBC(UBAR2, m_xnCells, m_ynCells, m_d);
    }
    else if(test_num == 7)//Kelvin test case
    {
      KelvintestBC(UBAR2, m_xnCells, m_ynCells, m_d, m_nVar);
    }
    else
    {
      TransBC(UBAR2, m_xnCells, m_ynCells, m_d);
    }

    // source term update
    for(int i = 2; i<m_xnCells+2; i++)
    {
      for(int j = 2; j<m_ynCells+2; j++)
      {
        for(int l = 0; l<m_nVar; l++)
        {
          U(i,j,0,l) =  UBAR2(i,j,0,l);
        }
          U(i,j,0,8) =  UBAR2(i,j,0,8)*exp(-dt*Ch*Ch/(0.18*Ch));
      }
    }
  }

  //convert u to prim
  for(int i=0; i<m_xnCells+4; i++)
  {
    for(int j=0; j<m_ynCells+4; j++)
    {
    E->set_gamma(m_gamma);
    PRIM(i,j,0) = E->u_to_prim(U(i,j,0),m_gamma,m_nVar,m_d);
    }
  }

  //outputting divergence cleaning data
  //L1 norm of the divergence nabla B 
  double L1 = 0;
  double nabla_B = 0;
  double nabla_B_max = 0;

  for(int i = 2; i<m_xnCells+2; i++)
  {
    for(int j = 2; j<m_ynCells+2; j++)
    {
      nabla_B = ((PRIM(i+1,j,0,7)-PRIM(i-1,j,0,7))/(2*m_dx)) + ((PRIM(i,j+1,0,5)-PRIM(i,j-1,0,5))/(2*m_dy));
      nabla_B_max = std::max(nabla_B_max, fabs(nabla_B));
      L1 = L1 + fabs(nabla_B);
    }
  }

  L1 = L1*dx*dy;

  std::cout<<L1<<std::endl;

  std::vector<double> test;
  test.resize(3);
  test[0] = t;
  test[1] = nabla_B_max;
  test[2] = L1;

  dc.push_back(test);

} while (t<m_tend);

//outputting nabla B 
std::ofstream output_DC ("DC.dat");
for(int k=0; k<counter; k++)
{
  output_DC<<dc[k][0]<<" "<<dc[k][1]<<" "<<dc[k][2]<<std::endl;
}
}

void system::outputData(int test_num)
{
//delete file from previous run
std::remove("MHD.dat"); 

//outputing data
std::ofstream output ("MHD.dat");
//std::ofstream dc ("conv.txt");
for(int i=2; i<m_xnCells+2; i++)
{
  for(int j=2; j<m_ynCells+2; j++)
  {
  double x = m_x0 + (i-1)*m_dx;
  double y = m_y0 + (j-1)*m_dy;
  output<<" "<<x<<" "<<y<<" "<<PRIM(i,j,0,0)<<" "<<PRIM(i,j,0,1)<<" "<<PRIM(i,j,0,2)<<" "<<PRIM(i,j,0,3)<<" "<<PRIM(i,j,0,4)<<" "<<PRIM(i,j,0,5)<<" "<<PRIM(i,j,0,6)<<" "<<PRIM(i,j,0,7)<<" "<<sqrt(pow(PRIM(i,j,0,7),2)+pow(PRIM(i,j,0,5),2))/PRIM(i,j,0,6)<<" "<<PRIM(i,j,0,4)/((m_gamma-1)*PRIM(i,j,0,0))<<std::endl;
  }
  output<<std::endl;
}
}
