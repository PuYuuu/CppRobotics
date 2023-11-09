#pragma once
#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <chrono>
#include <vector>
#include <string>
#include <numeric>
#include <climits>
#include <algorithm>

#include <Eigen/Core>

using namespace Eigen;

namespace utils {

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
Matrix2d rotation_matrix2d(T theta)
{
    Matrix2d rotation;
    rotation << cos(theta), -sin(theta), 
            sin(theta), cos(theta);
    
    return rotation;
}

template <typename T>
double pi_2_pi(T theta)
{
    while (theta > M_PI) {
        theta -= 2.0 * M_PI;
    }
    while (theta < -M_PI) {
        theta += 2.0 * M_PI;
    }

    return theta;
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

template <typename T>
std::vector<T> diff(const std::vector<T>& vec)
{
    std::vector<T> ret;
    for (size_t idx = 1; idx < vec.size(); ++idx) {
        ret.push_back(vec[idx] - vec[idx - 1]);
    }
    
    return ret;
}

template <typename T>
std::vector<T> cumsum(std::vector<T> vec)
{
    std::vector<T> output;
    T tmp = 0;
    for(size_t idx = 0; idx < vec.size(); ++idx) {
        tmp += vec[idx];
        output.push_back(tmp);
    }

    return output;
}

template <typename T>
int search_index(std::vector<T> nums, T target)
{
    int left = 0, right = nums.size() - 1;
    while(left <= right){
        int mid = (right - left) / 2 + left;
        int num = nums[mid];
        if (num == target) {
            return mid;
        } else if (num > target) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }

    return -1;
}

template <typename T>
double variance(const std::vector<T>& data) 
{
    if (data.empty()) {
        return 0.0;
    }

    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();

    double variance = 0.0;
    for (const double& value : data) {
        variance += pow(value - mean, 2);
    }
    variance /= data.size();

    return variance;
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

class VehicleConfig
{
public:
    double RF;  // [m] distance from rear to vehicle front end of vehicle
    double RB;  // [m] distance from rear to vehicle back end of vehicle
    double W;   // [m] width of vehicle
    double WD;  // [m] distance between left-right wheels
    double WB;  // [m] Wheel base
    double TR;  // [m] Tyre radius
    double TW;  // [m] Tyre width
    double MAX_STEER;

    // Default parameters
    VehicleConfig(double scale = 1.0) : RF(3.3), RB(0.8), W(2.4), WB(2.5), TR(0.44), TW(0.7), MAX_STEER(0.65) {
        WD = 0.7 * W;
        
        RF *= scale;
        RB *= scale;
        W  *= scale;
        WD *= scale;
        WB *= scale;
        TR *= scale;
        TW *= scale;
    }
    ~VehicleConfig() {}
};

void draw_vehicle(Vector3d state, double steer, VehicleConfig c,
                bool draw_wheel = true, std::string color="-k");

}

#endif
