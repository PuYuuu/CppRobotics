#pragma once
#ifndef __CUBIC_SPLINE_HPP
#define __CUBIC_SPLINE_HPP

#include <vector>
#include <Eigen/Eigen>

using std::vector;
using namespace Eigen;

class CubicSpline
{
private:
    vector<double> x;
    vector<double> y;
    vector<double> a;
    vector<double> b;
    vector<double> c;
    vector<double> d;
    vector<double> h;
    int nx;
public:
    CubicSpline() {}
    CubicSpline(vector<double> _x, vector<double> _y);
    ~CubicSpline() {}
    MatrixXd calc_A(void);
    VectorXd calc_B(void);
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
    vector<double> s;

    CubicSpline2D() {}
    CubicSpline2D(vector<double> _x, vector<double> _y);
    ~CubicSpline2D() {}
    vector<double> calc_s(vector<double> _x, vector<double> _y);
    Vector2d calc_position(double _s);
    double calc_yaw(double _s);
    double calc_curvature(double s);
};

#endif
