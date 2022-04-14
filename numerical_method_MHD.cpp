//This file contains the functions for slope reconstruction, wave speed estimate, HLLC and FORCE solver for solving MHD system
#include "numerical_method_MHD.H"

eos* E1 = new eos();

//constructor
numerical_method_MHD::numerical_method_MHD()
{}

//function to find Minbee constant
Arrayofdouble numerical_method_MHD::Minbee(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  for(int i=0; i<4; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
  }

  Arrayofdouble Mb;
  for(int i=0; i<4; i++)
  {
    if(r[i]<=0)
    {
      Mb[i] = 0;
    }
    else if(r[i]>0 && r[i]<=1)
    {
      Mb[i] = r[i];
    }
    else
    {
      if((2/(1+r[i])>1))
      {
        Mb[i] = 1;
      }
      else
      {
        Mb[i] = 2/(1+r[i]);
      }
    }
  }
  return Mb;
}


//function to find Vanleer constant
Arrayofdouble numerical_method_MHD::Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  Arrayofdouble xi;
  Arrayofdouble Vl;

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
    xi[i] = 2/(1+r[i]);
  }

  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Vl[i] = 0;
    }
    else if(r[i]>0)
    {
      if(xi[i]>(2*r[i]/(1+r[i])))
      {
        Vl[i] = 2*r[i]/(1+r[i]);
      }
      else 
      {
        Vl[i] = xi[i];
      }
    }
  }
  return Vl;
}

//defining wavespeedestimate for HLL and HLLC
//p represents primR_nplushalf[i], q represents primL_nplushalf[i+1]
//Sl(minimum wave speed from minimum vx - maximum Cf), Sr(maximum wave speed from maximum vx + maximum Cf), S_star (intermediate wave speed)
void numerical_method_MHD::wavespeedestimate(Arrayofdouble p, Arrayofdouble q, double n, double gamma, Arrayofdouble& wavespeed, int D)
{
  double rho_l = p[0];
  double rho_r = q[0];
  double vx_l = p[1];
  double vx_r = q[1];
  double vy_l = p[2];
  double vy_r = q[2];
  double vz_l = p[3];
  double vz_r = q[3];
  double p_l = p[4];
  double p_r = q[4];

  double Bx_l = p[7];
  double Bx_r = q[7];
  double By_l = p[5];
  double By_r = q[5];
  double Bz_l = p[6];
  double Bz_r = q[6];

  // compute wave speed Cf left
  // computing Ca
  double Ca_2xl = (Bx_l*Bx_l + By_l*By_l + Bz_l*Bz_l) / (rho_l);

  //computing Cs using total pressure
  double Cs_2xl = gamma*p_l/rho_l;

  //computing Cf
  double Cf_xl = sqrt(0.5*((Cs_2xl+Ca_2xl)+sqrt(pow((Cs_2xl+Ca_2xl),2)-(4*Cs_2xl*pow(Bx_l,2)/rho_l))));

  //compute wave speed Cf right
  //computing Ca
  double Ca_2xr = (Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r) / (rho_r);

  //computing Cs using total pressure
  double Cs_2xr = gamma*p_r/rho_r;

  //computing Cf
  double Cf_xr = sqrt(0.5*((Cs_2xr+Ca_2xr)+sqrt(pow((Cs_2xr+Ca_2xr),2)-(4*Cs_2xr*pow(Bx_r,2)/rho_r))));

  //compute Sl, Sr
  double Sl = std::min(vx_l, vx_r)-std::max(Cf_xl, Cf_xr);
  double Sr = std::max(vx_l, vx_r)+std::max(Cf_xl, Cf_xr);
  
  double Bx2_l = Bx_l*Bx_l;
  double Bx2_r = Bx_r*Bx_r;
  double B2_l = Bx_l*Bx_l + By_l*By_l + Bz_l*Bz_l;
  double B2_r = Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r;
  double S_star = (rho_r*vx_r*(Sr-vx_r) - rho_l*vx_l*(Sl-vx_l) + (p_l+0.5*B2_l) - (p_r+0.5*B2_r) - Bx2_l + Bx2_r) / (rho_r*(Sr-vx_r)-rho_l*(Sl-vx_l));

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;
}

//y wavespeed estimate
void numerical_method_MHD::ywavespeedestimate(Arrayofdouble p, Arrayofdouble q, double n, double gamma, Arrayofdouble& ywavespeed, int D)
{
  double rho_l = p[0];
  double rho_r = q[0];
  double vx_l = p[1];
  double vx_r = q[1];
  double vy_l = p[2];
  double vy_r = q[2];
  double vz_l = p[3];
  double vz_r = q[3];
  double p_l = p[4];
  double p_r = q[4];

  double Bx_l = p[7];
  double Bx_r = q[7];
  double By_l = p[5];
  double By_r = q[5];
  double Bz_l = p[6];
  double Bz_r = q[6];

  //compute wave speed Cf left
  //computing Ca
  double Ca_2yl = (Bx_l*Bx_l + By_l*By_l + Bz_l*Bz_l) / (rho_l);

  //computing Cs using total pressure
  double Cs_2yl = gamma*p_l/rho_l;

  //computing Cf
  double Cf_yl = sqrt(0.5*((Cs_2yl+Ca_2yl) + sqrt(pow((Cs_2yl+Ca_2yl),2)-(4*Cs_2yl*pow(By_l,2)/rho_l))));

  //compute wave speed Cf right
  //computing Ca
  double Ca_2yr = (Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r) / (rho_r);

  //computing Cs using total pressure
  double Cs_2yr = gamma*p_r/rho_r;

  //computing Cf
  double Cf_yr = sqrt(0.5*((Cs_2yr+Ca_2yr) + sqrt(pow((Cs_2yr+Ca_2yr),2)-(4*Cs_2yr*pow(By_r,2)/rho_r))));

  //computing velocity magnitude
  double v_ml = sqrt(vx_l*vx_l + vy_l*vy_l + vz_l*vz_l);
  double v_mr = sqrt(vx_r*vx_r + vy_r*vy_r + vz_r*vz_r);

  //compute Sl, Sr
  double Sl = std::min(vy_l, vy_r)-std::max(Cf_yl, Cf_yr);
  double Sr = std::max(vy_l, vy_r)+std::max(Cf_yl, Cf_yr);
  double By2_l = By_l*By_l;
  double By2_r = By_r*By_r;
  double B2_l = Bx_l*Bx_l + By_l*By_l + Bz_l*Bz_l;
  double B2_r = Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r;
  double S_star = (rho_r*vy_r*(Sr-vy_r) - rho_l*vy_l*(Sl-vy_l) + (p_l+0.5*B2_l) - (p_r+0.5*B2_r) - By2_l + By2_r) / (rho_r*(Sr-vy_r)-rho_l*(Sl-vy_l));

  ywavespeed[0] = Sl;
  ywavespeed[1] = Sr;
  ywavespeed[2] = S_star;
}

//computing nplushalf
void numerical_method_MHD::compute_nplushalf_variables(Arrayofdouble a, Arrayofdouble b, Arrayofdouble c, double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D)
{
  Arrayofdouble ui, uiMinus1, uiPlus1, Vl, Mb, delta;
  Arrayofdouble uiL, uiR, primiL, primiR, fuiL, fuiR;

  //defining ui, uiMinus1, uiPlus1
  for(int i=0; i<nVar; i++)
  {
  uiMinus1[i] = a[i];
  ui[i] = b[i];
  uiPlus1[i] = c[i];
  }

  //finding Vl
  Vl = Vanleer(ui, uiMinus1, uiPlus1, nVar);
  double Vl_min;
  Vl_min = *std::min_element(Vl.begin(), Vl.end());

  //finding Mb
  Mb = Minbee(ui, uiMinus1, uiPlus1, nVar);
  double Mb_min;
  Mb_min = *std::min_element(Mb.begin(), Mb.end());

  //define delta
  for(int i=0; i<nVar; i++)
  {
    delta[i] = 0.5*(uiPlus1[i]-uiMinus1[i]);
  }

  //define uiL, uiR
  for(int i=0; i<nVar; i++)
  {
    uiL[i] = ui[i] - 0.5*Mb_min*delta[i];
    uiR[i] = ui[i] + 0.5*Mb_min*delta[i];
    // uiL[i] = ui[i] - 0.5*Vl_min*delta[i];
    // uiR[i] = ui[i] + 0.5*Vl_min*delta[i];
  }

  //convert uiL, uiR to primiL, primiR
  E1->set_gamma(gamma);
  primiL = E1->u_to_prim(uiL, gamma, nVar, D);
  primiR = E1->u_to_prim(uiR, gamma, nVar, D);

  fuiL = E1->fluxf(uiL, primiL, nVar, D);
  fuiR = E1->fluxf(uiR, primiL, nVar, D);

  //calculate uiL_nplushalf, uiR_nplushalf
  for(int i=0; i<nVar; i++)
  {
    uiL_nplushalf[i] = uiL[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
    uiR_nplushalf[i] = uiR[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
  }
}

//computing nplushalf
void numerical_method_MHD::ycompute_nplushalf_variables(Arrayofdouble a, Arrayofdouble b, Arrayofdouble c, double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D)
{
  Arrayofdouble ui, uiMinus1, uiPlus1, Vl, Mb, delta;
  Arrayofdouble uiL, uiR, primiL, primiR, fuiL, fuiR;

  //defining ui, uiMinus1, uiPlus1
  for(int i=0; i<nVar; i++)
  {
  uiMinus1[i] = a[i];
  ui[i] = b[i];
  uiPlus1[i] = c[i];
  }

  //finding Vl
  Vl = Vanleer(ui, uiMinus1, uiPlus1, nVar);
  double Vl_min;
  Vl_min = *std::min_element(Vl.begin(), Vl.end());

  //finding Mb
  Mb = Minbee(ui, uiMinus1, uiPlus1, nVar);
  double Mb_min;
  Mb_min = *std::min_element(Mb.begin(), Mb.end());

  //define delta
  for(int i=0; i<nVar; i++)
  {
    delta[i] = 0.5*(uiPlus1[i]-uiMinus1[i]);
  }

  //define uiL, uiR
  for(int i=0; i<nVar; i++)
  {
    uiL[i] = ui[i] - 0.5*Mb_min*delta[i];
    uiR[i] = ui[i] + 0.5*Mb_min*delta[i];
    // uiL[i] = ui[i] - 0.5*Vl_min*delta[i];
    // uiR[i] = ui[i] + 0.5*Vl_min*delta[i];
  }

  //convert uiL, uiR to primiL, primiR
  E1->set_gamma(gamma);
  primiL = E1->u_to_prim(uiL, gamma, nVar, D);
  primiR = E1->u_to_prim(uiR, gamma, nVar, D);

  fuiL = E1->y_fluxf(uiL, primiL, nVar, D);
  fuiR = E1->y_fluxf(uiR, primiL, nVar, D);

  //calculate uiL_nplushalf, uiR_nplushalf
  for(int i=0; i<nVar; i++)
  {
    uiL_nplushalf[i] = uiL[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
    uiR_nplushalf[i] = uiR[i] - 0.5*dt/dx*(fuiR[i]-fuiL[i]);
  }
}

//divergence cleaning Bx term
double numerical_method_MHD::compute_Bx_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch)
{
  double Bx_l = a[7];
  double Bx_r = b[7];
  double psi_l = a[8];
  double psi_r = b[8];

  double Bx_dc = 0.5*(Bx_l + Bx_r) - (1/(2*Ch)) * (psi_r - psi_l);

  return Bx_dc;
}

//divergence cleaning By term
double numerical_method_MHD::compute_By_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch)
{
  double By_l = a[5];
  double By_r = b[5];
  double psi_l = a[8];
  double psi_r = b[8];

  double By_dc = 0.5*(By_l + By_r) - (1/(2*Ch)) * (psi_r - psi_l);

  return By_dc;
}

//divergence cleaning psi_x term
double numerical_method_MHD::compute_psix_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch)
{
  double Bx_l = a[7];
  double Bx_r = b[7];
  double psi_l = a[8];
  double psi_r = b[8];

  double psix_dc = 0.5*(psi_l + psi_r) - (Ch/2) * (Bx_r - Bx_l);

  return psix_dc;
}

//divergence cleaning By term
double numerical_method_MHD::compute_psiy_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch)
{
  double By_l = a[5];
  double By_r = b[5];
  double psi_l = a[8];
  double psi_r = b[8];

  double psiy_dc = 0.5*(psi_l + psi_r) - (Ch/2) * (By_r - By_l);

  return psiy_dc;
}

//divergence cleaning xflux
Arrayofdouble numerical_method_MHD::compute_xflux_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch)
{
  double Bx_l = a[7];
  double Bx_r = b[7];
  double psi_l = a[8];
  double psi_r = b[8];

  double Bx_dc = 0.5*(Bx_l + Bx_r) - (1/(2*Ch)) * (psi_r - psi_l);
  double psi = 0.5*(psi_l + psi_r) - (Ch/2) * (Bx_r - Bx_l);

  Arrayofdouble xflux_dc;

  xflux_dc[1] = Ch*Ch*Bx_dc;
  xflux_dc[0] = psi;

  return xflux_dc;
}

//divergence cleaning yflux
Arrayofdouble numerical_method_MHD::compute_yflux_dc(const Arrayofdouble& a, const Arrayofdouble& b, double Ch)
{
  double By_l = a[5];
  double By_r = b[5];
  double psi_l = a[8];
  double psi_r = b[8];

  double By_dc = 0.5*(By_l + By_r) - (1/(2*Ch)) * (psi_r - psi_l);
  double psi = 0.5*(psi_l + psi_r) - (Ch/2) * (By_r - By_l);

  Arrayofdouble yflux_dc;

  yflux_dc[1] = Ch*Ch*By_dc;
  yflux_dc[0] = psi;

  return yflux_dc;
}

//defining getHLLvar
void numerical_method_MHD::getHLLvar(Arrayofdouble a, Arrayofdouble b, Arrayofdouble p, Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& u_HLL, int nVar, int D, int div_clean)
{
  //calculate u_HLL to get Bx, By, Bz, vx, vy, vz for MHD HLLC solver
  Arrayofdouble flux_l; //fl
  Arrayofdouble flux_r; //fr

  double Bx_tilde = a[7];
  double psix_tilde = a[8];

  flux_l = E1->fluxf(a, p, nVar, D);
  flux_r = E1->fluxf(b, q, nVar, D);

  for(int i=0; i<nVar; i++)
  {
    u_HLL[i] = (wavespeed[1]*b[i]-wavespeed[0]*a[i]-(flux_r[i]-flux_l[i]))/(wavespeed[1]-wavespeed[0]);
  }

  if(div_clean == 1)
  {
    u_HLL[7] = Bx_tilde;//Bx input to HLLC needs to be BX tilde
    u_HLL[8] = psix_tilde;
  }
}

//calculate u_HLLC left (uL, uR, uL_star, uR_star) using fluxfunction and relationships derived from RH conditions
void numerical_method_MHD::getHLLCvar_l(Arrayofdouble& u_HLL, Arrayofdouble a, Arrayofdouble p, Arrayofdouble wavespeed, Arrayofdouble& HLLCvar_l, Arrayofdouble b, Arrayofdouble q, int nVar, int D)
{
  //define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  double rho_l = p[0];
  double vx_l = p[1];
  double vy_l = p[2];
  double vz_l = p[3];
  double p_l = p[4];
  double By_l = p[5];
  double Bz_l = p[6];
  double Bx_l = p[7];
  double U_l = a[4];

  double vx_star_l = S_star; 
  double vy_star_l = u_HLL[2] / u_HLL[0];
  double vz_star_l = u_HLL[3] / u_HLL[0];

  //compute u_Hllc
  double Bx_star_l = u_HLL[7];
  double By_star_l = u_HLL[5];
  double Bz_star_l = u_HLL[6];
  double rho_star_l = rho_l * (Sl - vx_l) / (Sl - S_star);
  double xmomentum_star_l = rho_l * (Sl - vx_l) / (Sl - S_star) * S_star;
  double ymomentum_star_l = rho_l * vy_l * (Sl - vx_l) / (Sl - S_star) - (Bx_star_l * By_star_l - Bx_l * By_l) / (Sl - S_star);
  double zmomentum_star_l = rho_l * vz_l * (Sl - vx_l) / (Sl - S_star) - (Bx_star_l * Bz_star_l - Bx_l * Bz_l) / (Sl - S_star);
  
  double B2_l = Bx_l*Bx_l + By_l*By_l + Bz_l*Bz_l;
  double p_star = rho_l * (Sl - vx_l) * (S_star - vx_l) + p_l + (0.5*B2_l) - (Bx_l * Bx_l) + (Bx_star_l * Bx_star_l); 

  //compute energy Hllc
  double Bv_star_l = Bx_star_l * (u_HLL[1] / u_HLL[0]) + By_star_l * vy_star_l + Bz_star_l * vz_star_l; 
  double Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;

  double U_star_l = U_l * (Sl - vx_l) / (Sl - S_star) + (p_star * S_star - (p_l+0.5*B2_l) * vx_l - (Bx_star_l * Bv_star_l - Bx_l * Bv_l)) / (Sl - vx_star_l);

  //define an array of HLLC variables
  HLLCvar_l[0] = rho_star_l;
  HLLCvar_l[1] = xmomentum_star_l;
  HLLCvar_l[2] = ymomentum_star_l;
  HLLCvar_l[3] = zmomentum_star_l;
  HLLCvar_l[4] = U_star_l;
  HLLCvar_l[5] = By_star_l;
  HLLCvar_l[6] = Bz_star_l;
  HLLCvar_l[7] = Bx_star_l;
}

//calculate u_HLLC right (uL, uR, uL_star, uR_star) using fluxfunction and relationships derived from RH conditions
void numerical_method_MHD::getHLLCvar_r(Arrayofdouble& u_HLL, Arrayofdouble b, Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& HLLCvar_r, Arrayofdouble a, Arrayofdouble p, int nVar, int D)
{
  //define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  double rho_r = q[0];
  double vx_r = q[1];
  double vy_r = q[2];
  double vz_r = q[3];
  double p_r = q[4];
  double By_r = q[5];
  double Bz_r = q[6];
  double Bx_r = q[7];
  double U_r = b[4];

  double vx_star_r = S_star; 
  double vy_star_r = u_HLL[2] / u_HLL[0];
  double vz_star_r = u_HLL[3] / u_HLL[0];

  //compute u_Hllc
  double Bx_star_r = u_HLL[7]; //Bx_Hllc equal to Bx_l, Bx_r, Bx_Hll in 1-D
  double By_star_r = u_HLL[5];
  double Bz_star_r = u_HLL[6];
  double rho_star_r = rho_r * (Sr - vx_r) / (Sr - S_star);
  double xmomentum_star_r = rho_star_r * S_star;
  double ymomentum_star_r = rho_r * vy_r * (Sr - vx_r) / (Sr - S_star) - (Bx_star_r * By_star_r - Bx_r * By_r) / (Sr - S_star);
  double zmomentum_star_r = rho_r * vz_r * (Sr - vx_r) / (Sr - S_star) - (Bx_star_r * Bz_star_r - Bx_r * Bz_r) / (Sr - S_star);
  double B2_r = Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r;
  double p_star = rho_r * (Sr - vx_r) * (S_star - vx_r) + p_r + (0.5*B2_r) - (Bx_r * Bx_r) + (Bx_star_r * Bx_star_r); //p_star_l equals to p_star_r

  //compute energy Hllc
  double Bv_star_r = Bx_star_r * (u_HLL[1] / u_HLL[0])  + By_star_r * vy_star_r + Bz_star_r * vz_star_r;
  double Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
  double U_star_r = U_r * (Sr - vx_r) / (Sr - S_star) + (p_star * S_star - (p_r+0.5*B2_r) * vx_r - (Bx_star_r * Bv_star_r - Bx_r * Bv_r)) / (Sr - vx_star_r);

  //define an array of HLLC variables
  HLLCvar_r[0] = rho_star_r;
  HLLCvar_r[1] = xmomentum_star_r;
  HLLCvar_r[2] = ymomentum_star_r;
  HLLCvar_r[3] = zmomentum_star_r;
  HLLCvar_r[4] = U_star_r;
  HLLCvar_r[5] = By_star_r;
  HLLCvar_r[6] = Bz_star_r;
  HLLCvar_r[7] = Bx_star_r;
}

//defining ygetHLLvar
void numerical_method_MHD::ygetHLLvar(Arrayofdouble a, Arrayofdouble b, Arrayofdouble p, Arrayofdouble q, Arrayofdouble ywavespeed, Arrayofdouble& yu_HLL, int nVar, int D, int div_clean)
{
  //calculate u_HLL to get Bx, By, Bz, vx, vy, vz for MHD HLLC solver
  Arrayofdouble flux_l; 
  Arrayofdouble flux_r; 

  double By_tilde = a[5];
  double psiy_tilde = a[8];

  flux_l = E1->y_fluxf(a, p, nVar, D);
  flux_r = E1->y_fluxf(b, q, nVar, D);

  for(int i=0; i<nVar; i++)
  {
    yu_HLL[i] = (ywavespeed[1]*b[i]-ywavespeed[0]*a[i]-(flux_r[i]-flux_l[i]))/(ywavespeed[1]-ywavespeed[0]);
  }

  if(div_clean == 1)
  {
    yu_HLL[5] = By_tilde;//By input into HLLC solver should be By tilde
    yu_HLL[8] = psiy_tilde;
  }
}

//HLLC variables for y direction
//calculate u_HLLC left (uL, uR, uL_star, uR_star) using fluxfunction and relationships derived from RH conditions
void numerical_method_MHD::ygetHLLCvar_l(Arrayofdouble& yu_HLL, Arrayofdouble a, Arrayofdouble p, Arrayofdouble ywavespeed, Arrayofdouble& yHLLCvar_l, Arrayofdouble b, Arrayofdouble q, int nVar, int D)
{
  //define input variables 
  double Sl = ywavespeed[0];
  double Sr = ywavespeed[1];
  double S_star = ywavespeed[2];

  double rho_l = p[0];
  double vx_l = p[1];
  double vy_l = p[2];
  double vz_l = p[3];
  double p_l = p[4];
  double By_l = p[5];
  double Bz_l = p[6];
  double Bx_l = p[7];
  double U_l = a[4];

  double vx_star_l = S_star; 
  double vy_star_l = yu_HLL[2] / yu_HLL[0];
  double vz_star_l = yu_HLL[3] / yu_HLL[0];

  //compute u_Hllc
  double Bx_star_l = yu_HLL[7];
  double By_star_l = yu_HLL[5];
  double Bz_star_l = yu_HLL[6];
  double rho_star_l = rho_l * (Sl - vy_l) / (Sl - S_star);
  double xmomentum_star_l = rho_star_l * vx_l - (Bx_star_l * By_star_l - Bx_l * By_l) / (Sl - S_star);
  double ymomentum_star_l = rho_l * S_star * (Sl - vy_l) / (Sl - S_star);
  double zmomentum_star_l = rho_l * vz_l * (Sl - vy_l) / (Sl - S_star) - (By_star_l * Bz_star_l - By_l * Bz_l) / (Sl - S_star);

  double B2_l = Bx_l*Bx_l + By_l*By_l + Bz_l*Bz_l;
  double p_star = rho_l * (Sl - vy_l) * (S_star - vy_l) + p_l + (0.5*B2_l)  - (By_l * By_l) + (By_star_l * By_star_l);

  //compute energy Hllc
  double Bv_star_l = Bx_star_l * (yu_HLL[1] / yu_HLL[0])+ By_star_l * vy_star_l + Bz_star_l * vz_star_l;
  double Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;
  double U_star_l = U_l * (Sl - vy_l) / (Sl - S_star) + (p_star * S_star - (p_l+0.5*B2_l) * vy_l - (By_star_l * Bv_star_l - By_l * Bv_l)) / (Sl - S_star);

  //define an array of HLLC variables
  yHLLCvar_l[0] = rho_star_l;
  yHLLCvar_l[1] = xmomentum_star_l;
  yHLLCvar_l[2] = ymomentum_star_l;
  yHLLCvar_l[3] = zmomentum_star_l;
  yHLLCvar_l[4] = U_star_l;
  yHLLCvar_l[5] = By_star_l;
  yHLLCvar_l[6] = Bz_star_l;
  yHLLCvar_l[7] = Bx_star_l;
}

//calculate u_HLLC right (uL, uR, uL_star, uR_star) using fluxfunction and relationships derived from RH conditions
void numerical_method_MHD::ygetHLLCvar_r(Arrayofdouble& yu_HLL, Arrayofdouble b, Arrayofdouble q, Arrayofdouble ywavespeed, Arrayofdouble& yHLLCvar_r, Arrayofdouble a, Arrayofdouble p, int nVar, int D)
{
  //define input variables 
  double Sl = ywavespeed[0];
  double Sr = ywavespeed[1];
  double S_star = ywavespeed[2];

  double rho_r = q[0];
  double vx_r = q[1];
  double vy_r = q[2];
  double vz_r = q[3];
  double p_r = q[4];
  double By_r = q[5];
  double Bz_r = q[6];
  double Bx_r = q[7];
  double U_r = b[4];

  double vx_star_r = S_star;
  double vy_star_r = yu_HLL[2] / yu_HLL[0];
  double vz_star_r = yu_HLL[3] / yu_HLL[0];

  //compute u_Hllc
  double Bx_star_r = yu_HLL[7]; //Bx_Hllc equal to Bx_l, Bx_r, Bx_Hll in 1-D
  double By_star_r = yu_HLL[5];
  double Bz_star_r = yu_HLL[6];
  double rho_star_r = rho_r * (Sr - vy_r) / (Sr - S_star);

  double xmomentum_star_r = rho_star_r * vx_r - (Bx_star_r * By_star_r - Bx_r * By_r) / (Sr - S_star);
  double ymomentum_star_r = rho_r * S_star * (Sr - vy_r) / (Sr - S_star);
  double zmomentum_star_r = rho_r * vz_r * (Sr - vy_r) / (Sr - S_star) - (By_star_r * Bz_star_r - By_r * Bz_r) / (Sr - S_star);

  double B2_r = Bx_r*Bx_r + By_r*By_r + Bz_r*Bz_r;
  double p_star = rho_r * (Sr - vy_r) * (S_star - vy_r) + p_r + (0.5*B2_r) - (By_r * By_r) + (By_star_r * By_star_r);

  //compute energy Hllc
  double Bv_star_r = Bx_star_r * (yu_HLL[1] / yu_HLL[0]) + By_star_r * vy_star_r + Bz_star_r * vz_star_r;
  double Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
  double U_star_r = U_r * (Sr - vy_r) / (Sr - S_star) + (p_star * S_star - (p_r+0.5*B2_r) * vy_r - (By_star_r * Bv_star_r - By_r * Bv_r)) / (Sr - S_star);

  //define an array of HLLC variables
  yHLLCvar_r[0] = rho_star_r;
  yHLLCvar_r[1] = xmomentum_star_r;
  yHLLCvar_r[2] = ymomentum_star_r;
  yHLLCvar_r[3] = zmomentum_star_r;
  yHLLCvar_r[4] = U_star_r;
  yHLLCvar_r[5] = By_star_r;
  yHLLCvar_r[6] = Bz_star_r;
  yHLLCvar_r[7] = Bx_star_r;
}

//defining a function for HLLC Xflux
void numerical_method_MHD::compute_HLLCflux(Arrayofdouble a, Arrayofdouble p, Arrayofdouble b, Arrayofdouble q, Arrayofdouble wavespeed, Arrayofdouble& Fhllc, int nVar, double gamma, int D, Arrayofdouble& u_HLL, Arrayofdouble& HLLCvar_l, Arrayofdouble& HLLCvar_r, int div_clean, double Ch)
{
  Arrayofdouble xflux_dc;
  if(div_clean == 1)
  {
  //compute Bx_dc for divergence cleaning
  double Bx_dc = compute_Bx_dc(a,b,Ch);

  //compute psix_dc for divergence cleaning
  double psix_dc = compute_psix_dc(a,b,Ch);

  //compute xflux for divergence cleaning
  xflux_dc = compute_xflux_dc(a,b,Ch);

  //modifying data to be passed into computing HLL variables
  a[7] = Bx_dc;
  b[7] = Bx_dc;
  p[7] = Bx_dc;
  q[7] = Bx_dc;

  a[8] = psix_dc;
  b[8] = psix_dc;
  p[8] = psix_dc;
  q[8] = psix_dc;
  }

  //get HLL variables
  getHLLvar(a,b,p,q,wavespeed,u_HLL,nVar,D,div_clean);

  //get HLLC left variables
  getHLLCvar_l(u_HLL, a, p, wavespeed, HLLCvar_l, b, q, nVar, D);

  //get HLLC right variables
  getHLLCvar_r(u_HLL, b, q, wavespeed, HLLCvar_r, a, p, nVar, D);

  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];
  Arrayofdouble HLLCprim_l, HLLCprim_r, Fl, Fr, Fl_star, Fr_star;

  HLLCprim_l = E1->u_to_prim(HLLCvar_l, gamma, nVar, D);
  HLLCprim_r = E1->u_to_prim(HLLCvar_r, gamma, nVar, D);
  Fl = E1->fluxf(a, p, nVar, D);
  Fr = E1->fluxf(b, q, nVar, D);

  for(int i=0; i<nVar; i++)
  {
  Fl_star[i] = Fl[i]  + Sl*(HLLCvar_l[i] - a[i]);
  Fr_star[i] = Fr[i]  + Sr*(HLLCvar_r[i] - b[i]);
  }

  //logic for choosing the correct flux
  if(div_clean == 0)
  {
    if(Sl >= 0)
    {	
        Fhllc = Fl;
        return; 
    }

    if(Sl < 0 && S_star >= 0)
    {
        Fhllc = Fl_star;
        return;
    }

    if(S_star < 0 && Sr > 0)
    {
        Fhllc = Fr_star;
        return;
    }

    if(Sr <= 0)
    {	
        Fhllc = Fr;
        return; 
    }	
  }
  else if(div_clean == 1)	
  {
     //flux for Bx and psi can be computed separately
    if(Sl >= 0)
    {	
        Fhllc = Fl;
        Fhllc[7] = xflux_dc[0];
        Fhllc[8] = xflux_dc[1];
        return; 
    }

    if(Sl < 0 && S_star >= 0)
    {
        Fhllc = Fl_star;
        Fhllc[7] = xflux_dc[0];
        Fhllc[8] = xflux_dc[1];
        return;
    }

    if(S_star < 0 && Sr > 0)
    {
        Fhllc = Fr_star;
        Fhllc[7] = xflux_dc[0];
        Fhllc[8] = xflux_dc[1];
        return;
    }

    if(Sr <= 0)
    {	
        Fhllc = Fr;
        Fhllc[7] = xflux_dc[0];
        Fhllc[8] = xflux_dc[1];
        return; 
    }		 
  }
}

//defining a function to get HLLC Yflux
void numerical_method_MHD::ycompute_HLLCflux(Arrayofdouble a, Arrayofdouble p, Arrayofdouble b, Arrayofdouble q, Arrayofdouble ywavespeed, Arrayofdouble& Fhllc, int nVar, double gamma, int D, Arrayofdouble& yu_HLL, Arrayofdouble& yHLLCvar_l, Arrayofdouble& yHLLCvar_r, int div_clean, double Ch)
{
  Arrayofdouble yflux_dc;
  if(div_clean == 1)
  {
  //compute By_dc for divergence cleaning
  double By_dc = compute_By_dc(a,b,Ch);

  //compute psiy_dc for divergence cleaning
  double psiy_dc = compute_psiy_dc(a,b,Ch);

  //compute yflux for divergence cleaning
  yflux_dc = compute_yflux_dc(a,b,Ch);

  //modifying data to be passed into computing HLL variables
  a[5] = By_dc;
  b[5] = By_dc;
  p[5] = By_dc;
  q[5] = By_dc;

  a[8] = psiy_dc;
  b[8] = psiy_dc;
  p[8] = psiy_dc;
  q[8] = psiy_dc;
  }

  //get HLL variables
  ygetHLLvar(a,b,p,q,ywavespeed,yu_HLL,nVar,D,div_clean);

  //get HLLC left variables
  ygetHLLCvar_l(yu_HLL, a, p, ywavespeed, yHLLCvar_l, b, q, nVar, D);

  //get HLLC right variables
  ygetHLLCvar_r(yu_HLL, b, q, ywavespeed, yHLLCvar_r, a, p, nVar, D);

  double Sl = ywavespeed[0];
  double Sr = ywavespeed[1];
  double S_star = ywavespeed[2];
  Arrayofdouble HLLCprim_l, HLLCprim_r, Fl, Fr, Fl_star, Fr_star;

  HLLCprim_l = E1->u_to_prim(yHLLCvar_l, gamma, nVar, D);
  HLLCprim_r = E1->u_to_prim(yHLLCvar_r, gamma, nVar, D);
  Fl = E1->y_fluxf(a, p, nVar, D);
  Fr = E1->y_fluxf(b, q, nVar, D);

  for(int i=0; i<nVar; i++)
  {
  Fl_star[i] = Fl[i]  + Sl*(yHLLCvar_l[i] - a[i]);
  Fr_star[i] = Fr[i]  + Sr*(yHLLCvar_r[i] - b[i]);
  }

  if(div_clean == 0)
  {
    //logic for choosing the correct flux
    if(Sl >= 0)
    {	
        Fhllc = Fl;
        return; 
    }

    if(Sl < 0 && S_star >= 0)
    {
        Fhllc = Fl_star;
        return;
    }

    if(S_star < 0 && Sr > 0)
    {
        Fhllc = Fr_star;
        return;
    }

    if(Sr <= 0)
    {	Fhllc = Fr;
        return; 
    }		
  }
  else if(div_clean == 1)
  {
    if(Sl >= 0)
    {	
        Fhllc = Fl;
        Fhllc[5] = yflux_dc[0];
        Fhllc[8] = yflux_dc[1];
        return; 
    }

    if(Sl < 0 && S_star >= 0)
    {
        Fhllc = Fl_star;
        Fhllc[5] = yflux_dc[0];
        Fhllc[8] = yflux_dc[1];
        return;
    }

    if(S_star < 0 && Sr > 0)
    {
        Fhllc = Fr_star;
        Fhllc[5] = yflux_dc[0];
        Fhllc[8] = yflux_dc[1];
        return;
    }

    if(Sr <= 0)
    {	
        Fhllc = Fr;
        Fhllc[5] = yflux_dc[0];
        Fhllc[8] = yflux_dc[1];
        return; 
    }		
  }
}

//defining getFORCEflux for single cell (q=prim; p=prim(i+1))
Arrayofdouble numerical_method_MHD::getFORCEflux(Arrayofdouble a, Arrayofdouble b, Arrayofdouble p, Arrayofdouble q, double dx, double dt, double gamma, double nVar, int D)
{
  Arrayofdouble fluxfunc;
  Arrayofdouble fluxfuncplus1;
  Arrayofdouble LFflux;
  Arrayofdouble RIflux;
  Arrayofdouble RI;
  Arrayofdouble FORCEflux;

  E1->set_gamma(gamma);
  fluxfunc = E1->fluxf(a, p, nVar, D);
  fluxfuncplus1 = E1->fluxf(b, q, nVar, D);

  //FORCE method
  for(int k=0; k<nVar; k++)
  {
  LFflux[k] = (0.5*(dx/dt)*(a[k]-b[k]))+ (0.5*(fluxfunc[k]+fluxfuncplus1[k]));
  }

  // //RI in u(i+0.5)
  for(int i=0; i<nVar; i++)
  {
  RI[i] = (0.5*(a[i]+b[i])) - (0.5*(dt/dx)*(fluxfuncplus1[i]-fluxfunc[i]));
  }

  //converting RI using the same operations of u(i+0.5) to prim(i+0.5) conversion
  Arrayofdouble RI_prim;
  RI_prim = E1->u_to_prim(RI, gamma, nVar, D);

  //RI*prim(i+0.5)
  RIflux = E1->fluxf(RI, RI_prim, nVar, D);

  for(int i=0; i<nVar; i++)
  {
  FORCEflux[i] = 0.5*(LFflux[i]+RIflux[i]);
  }
  return FORCEflux;
}

//defining getFORCEflux for single cell (q=prim; p=prim(i+1))
Arrayofdouble numerical_method_MHD::ygetFORCEflux(Arrayofdouble a, Arrayofdouble b,Arrayofdouble p,Arrayofdouble q, double dx, double dt, double gamma, double nVar, int D)
{
  Arrayofdouble fluxfunc;
  Arrayofdouble fluxfuncplus1;
  Arrayofdouble LFflux;
  Arrayofdouble RIflux;
  Arrayofdouble RI;
  Arrayofdouble FORCEflux;

  E1->set_gamma(gamma);
  fluxfunc = E1->y_fluxf(a, p, nVar, D);
  fluxfuncplus1 = E1->y_fluxf(b, q, nVar, D);

  //FORCE method
  for(int k=0; k<nVar; k++)
  {
  LFflux[k] = (0.5*(dx/dt)*(a[k]-b[k]))+ (0.5*(fluxfunc[k]+fluxfuncplus1[k]));
  }

  // //RI in u(i+0.5)
  for(int i=0; i<nVar; i++)
  {
  RI[i] = (0.5*(a[i]+b[i])) - (0.5*(dt/dx)*(fluxfuncplus1[i]-fluxfunc[i]));
  }

  //converting RI using the same operations of u(i+0.5) to prim(i+0.5) conversion
  Arrayofdouble RI_prim;
  RI_prim = E1->u_to_prim(RI, gamma, nVar, D);

  //RI*prim(i+0.5)
  RIflux = E1->y_fluxf(RI, RI_prim, nVar, D);

  for(int i=0; i<nVar; i++)
  {
  FORCEflux[i] = 0.5*(LFflux[i]+RIflux[i]);
  }
  return FORCEflux;
}