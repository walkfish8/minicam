#include <minicam.hpp>

int main()
{
    minicam::Rotation_<float> r(1, 2, 1);
    std::cout << r << std::endl;
    return 0;
}