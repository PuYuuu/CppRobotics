#pragma once

#include <Eigen/Core>

using namespace Eigen;

namespace Utils {

template <typename T>
int sign(T num) {
    if (num < 0) {
        return -1;
    }

    return 1;
}

template <typename T>
Matrix3d transformation_matrix2d(T x, T y, T theta)
{
    Matrix3d trans;
    trans << cos(theta), -sin(theta), x, 
            sin(theta), cos(theta), y, 
            0, 0, 1;
    
    return trans;
}

}