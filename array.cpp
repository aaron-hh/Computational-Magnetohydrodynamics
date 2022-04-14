//This file contains the functions for creating the data structure required in the MHD simulations
#include "array.H"

m_array::m_array()
{}

void m_array::setSize(double XSIZE, double YSIZE, double ZSIZE)
{
    m_xSize = XSIZE;
    m_ySize = YSIZE;
    m_zSize = ZSIZE;
    m_TSize = XSIZE * YSIZE * ZSIZE;
    m_data.resize(m_TSize);
}

std::array<double,10>& m_array::operator()(int i, int j, int k) 
{
    return m_data[i + j*m_xSize + k*m_xSize*m_ySize];
}

double& m_array::operator()(int i, int j, int k, int l) 
{
    return m_data[i + j*m_xSize + k*m_xSize*m_ySize][l];
}
