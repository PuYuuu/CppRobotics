#pragma once

#include <chrono>
#include <vector>
#include <climits>
#include <algorithm>

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

template <typename T>
T max(std::vector<T> vec)
{
    int size = vec.size();
    assert(size > 0);

    T ret = vec[0];
    for (int idx = 1; idx < size; ++idx) {
        if (vec[idx] > ret) {
            ret = vec[idx];
        }
    }
    
    return ret;
}

template <typename T>
T min(std::vector<T> vec)
{
    int size = vec.size();
    assert(size > 0);

    T ret = vec[0];
    for (int idx = 1; idx < size; ++idx) {
        if (vec[idx] < ret) {
            ret = vec[idx];
        }
    }
    
    return ret;
}

class TicToc
{
public:
    TicToc() {
        tic();
    }

    void tic() {
        start = std::chrono::system_clock::now();
    }

    double toc() {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() * 1000;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> start, end;
};

}