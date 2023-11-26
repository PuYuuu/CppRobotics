#pragma once
#ifndef __CUBIC_SPLINE_HPP
#define __CUBIC_SPLINE_HPP

#include <vector>
#include <Eigen/Eigen>

class CubicSpline
{
private:
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
    std::vector<double> h;
    int nx;
public:
    CubicSpline() {}
    CubicSpline(std::vector<double> _x, std::vector<double> _y);
    ~CubicSpline() {}
    Eigen::MatrixXd calc_A(void);
    Eigen::VectorXd calc_B(void);
    double calc_position(double _x);
    double calc_first_derivative(double _x);
    double calc_second_derivative(double _x);
};

class CubicSpline2D
{
private:
    CubicSpline sx;
    CubicSpline sy;
public:
    std::vector<double> s;

    CubicSpline2D() {}
    CubicSpline2D(std::vector<double> _x, std::vector<double> _y);
    ~CubicSpline2D() {}
    std::vector<double> calc_s(std::vector<double> _x, std::vector<double> _y);
    Eigen::Vector2d calc_position(double _s);
    double calc_yaw(double _s);
    double calc_curvature(double s);
};

#endif
