#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <array>

class m_array
{
private:
    size_t m_xSize;
    size_t m_ySize;
    size_t m_zSize;
    size_t m_TSize;

    std::vector<std::array<double,10> >m_data; 

public:
    m_array();

    void setSize(double XSIZE, double YSIZE, double ZSIZE);
    std::array<double,10>& operator()(int i, int j, int k);
    double& operator()(int i, int j, int k, int l);
};


