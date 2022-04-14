//This file contains the function to calculate error and convergence rate

//read in files (100, 200, 400, 800, 8000) - run simulations manually at different resolutions
//calculate L1, L2, Linf norms for every resolutions
//calculate L1, L2, Linf orders for every resolutions - L1 order will be zero as previous error needed

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

typedef std::vector<double> Vectorofdouble;
typedef std::vector<std::vector<double> > Vector_vectorofdouble;

//defining global variables
int num_L = 3;

void compute_error(Vectorofdouble& Error, Vectorofdouble num_sol, Vectorofdouble exact_sol, double dx, int N)
{
    double sum1, sum2;
	
	double LInfnorm = 0;

	for(int i=0; i<N; i++)
	{
	sum1 += fabs(num_sol[i]-exact_sol[i*10*dx]);
	sum2 += pow((num_sol[i]-exact_sol[i*10*dx]),2);
	LInfnorm = std::max(LInfnorm,fabs(num_sol[i]-exact_sol[i*10*dx]));	
	}

	double L1norm = 0;
	L1norm = sum1/N;
	
	double L2norm = 0;
	L2norm += sqrt(sum2/N);

	Error[0] = L1norm;
	Error[1] = L2norm;
	Error[2] = LInfnorm;

}

void compute_convrate(Vector_vectorofdouble& ConvRates,Vector_vectorofdouble Error, int num_mesh)
{
    for(int i = 0; i<num_mesh-1; i++)//4 meshes in total
    {
        for(int j = 0; j<num_L; j++)//3 errors in total
        {
            ConvRates[0][j] = 0;
            ConvRates[i+1][j] = log(Error[i+1][j]/Error[i][j])/log(0.5);
        }
    }   
}

int main()
{
    //read in files containing numerical solution for 800 cells
    //std::ifstream ns_e("num_sol_eight.txt");
    std::ifstream ns_e("slic_eight.txt");

    Vectorofdouble v_e;

    double data1_e, data2_e;

    while(ns_e >> data1_e && ns_e >> data2_e)
    {
        v_e.push_back(data2_e);
    }

    ns_e.close();

    //read in files containing numerical solution for 400 cells
    //std::ifstream ns_f("num_sol_four.txt");
    std::ifstream ns_f("slic_four.txt");

    Vectorofdouble v_f;

    double data1_f, data2_f;

    while(ns_f >> data1_f && ns_f >> data2_f)
    {
        v_f.push_back(data2_f);
    }

    ns_f.close();

    //read in files containing numerical solution for 200 cells
    //std::ifstream ns_t("num_sol_two.txt");
    std::ifstream ns_t("slic_two.txt");

    Vectorofdouble v_t;

    double data1_t, data2_t;

    while(ns_t >> data1_t && ns_t >> data2_t)
    {
        v_t.push_back(data2_t);
    }

    ns_t.close();

    //read in files containing numerical solution for 100 cells
    //std::ifstream ns_o("num_sol_one.txt");
    std::ifstream ns_o("slic_one.txt");

    Vectorofdouble v_o;

    double data1_o, data2_o;

    while(ns_o >> data1_o && ns_o >> data2_o)
    {
        v_o.push_back(data2_o);
    }

    ns_o.close();

    //read in files containing exact solution 
    std::ifstream es("exact_sol.txt");

    Vectorofdouble v_exact;

    double data1_exact, data2_exact;

    while(es >> data1_exact && es >> data2_exact)
    {
        v_exact.push_back(data2_exact);
    }

    es.close();

    //calculate erros for all resolutions
    Vectorofdouble mesh(4);
    mesh[0] = 100;
    mesh[1] = 200;
    mesh[2] = 400;
    mesh[3] = 800;

    int num_mesh = mesh.size();

    Vector_vectorofdouble error_all;
    error_all.resize(num_mesh, Vectorofdouble(num_L));

    // 100 cells
    compute_error(error_all[0], v_o, v_exact, 8.0, 100);
    // 200 cells
    compute_error(error_all[1], v_t, v_exact, 4.0, 200);
    // 400 cells
    compute_error(error_all[2], v_f, v_exact, 2.0, 400);
    // 800 cells
    compute_error(error_all[3], v_e, v_exact, 1.0, 800);

    //calculate orders for all resolutions
    Vector_vectorofdouble order_all;
    order_all.resize(num_mesh, Vectorofdouble(num_L));
    compute_convrate(order_all, error_all, num_mesh);

    //output
    std::ofstream output_convergence ("CR.dat");

    output_convergence << "# Mesh - L1Error - L1ORDER - L2Error - L2ORDER - LInfError - LINFORDER" << std::endl;
            

    for(int i=0; i<num_mesh; i++)
    {
        output_convergence<<mesh[i]<<" "<<error_all[i][0]<<" "<<order_all[i][0]<<" "<<error_all[i][1]<<" "<<order_all[i][1]<<" "<<error_all[i][2]<<" "<<order_all[i][2]<<std::endl;
    }
}
