#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::vector;
namespace plt = matplotlibcpp;

unsigned long long factorial(int n)
{
    unsigned long long result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }

    return result;
}

unsigned long long comb(int n, int k) 
{
    if (k < 0 || k > n) {
        return 0;
    }

    return factorial(n) / (factorial(k) * factorial(n - k));
}

double bernstein_poly(int n, int i, double t)
{
    return comb(n, i) * pow(t, i) * pow(1 - t, n - i);
}

vector<double> bezier(double t, const vector<vector<double>>& control_points)
{
    int n = control_points[0].size() - 1;
    vector<double> point = {0.0, 0.0};
    for (size_t idx = 0; idx < (n + 1); ++idx) {
        double poly = bernstein_poly(n, idx, t);
        point[0] += (poly * control_points[0][idx]);
        point[1] += (poly * control_points[1][idx]);
    }

    return point;
}

vector<vector<double>> calc_bezier_path(const vector<vector<double>>& control_points,
    int n_points = 100)
{
    vector<vector<double>> traj(2);
    double step = 1.0 / n_points;
    for (double t = 0.0; t <= 1.0; t += step) {
        vector<double> point = bezier(t, control_points);
        traj[0].push_back(point[0]);
        traj[1].push_back(point[1]);
    }

    return traj;
}

vector<vector<double>> calc_4points_bezier_path(double sx, double sy, double syaw,
    double ex, double ey, double eyaw, double offset, vector<vector<double>>& cp)
{
    double dist = hypot(sx - ex, sy - ey) / offset;
    vector<vector<double>> control_points = {
        {sx, sx + dist * cos(syaw), ex - dist * cos(eyaw), ex},
        {sy, sy + dist * sin(syaw), ey - dist * sin(eyaw), ey}};

    cp = control_points;
    vector<vector<double>> path = calc_bezier_path(control_points);
    
    return path;
}

void dynamic_effect(void)
{
    vector<double> sx = {-3, 0, 4, 6};
    vector<double> sy = {2, 0, 1.5, 6};
    vector<double> pathx;
    vector<double> pathy;
    
    for (double t = 0; t <= 1.0; t += 0.01) {
        vector<double> x;
        vector<double> y;
        for (size_t idx = 0; idx < (sx.size() - 1); ++idx) {
            x.push_back(sx[idx + 1] * t + sx[idx] * (1 - t));
            y.push_back(sy[idx + 1] * t + sy[idx] * (1 - t));
        }

        vector<double> xx;
        vector<double> yy;
        for (size_t idx = 0; idx < (x.size() - 1); ++idx) {
            xx.push_back(x[idx + 1] * t + x[idx] * (1 - t));
            yy.push_back(y[idx + 1] * t + y[idx] * (1 - t));
        }

        double px = xx[1] * t + xx[0] * (1 - t);
        double py = yy[1] * t + yy[0] * (1 - t);
        pathx.push_back(px);
        pathy.push_back(py);

        plt::cla();
        plt::named_plot("Control Points", sx, sy, "-o");
        plt::plot(x, y);
        plt::plot(xx, yy);
        plt::named_plot("Bezier Path", pathx, pathy);
        plt::plot({px}, {py}, "o");
        plt::axis("equal");
        plt::legend();
        plt::title("Cubic Bezier Curve demo");
        plt::grid(true);
        plt::pause(0.001);
    }

    plt::show();
}

void bezier_path(void)
{
    double start_x = 10.0;
    double start_y = 1.0;
    double start_yaw = M_PI;

    double end_x = -0.0;
    double end_y = -3.0;
    double end_yaw = -M_PI_4;
    double offset = 3.0;
    vector<vector<double>> control_points;
    vector<vector<double>> path = calc_4points_bezier_path(start_x, start_y,
                    start_yaw, end_x, end_y, end_yaw, offset, control_points);
    
    plt::named_plot("Bezier Path", path[0], path[1]);
    plt::named_plot("Control Points", control_points[0], control_points[1], "--o");
    plt::arrow(start_x, start_y, cos(start_yaw), sin(start_yaw), "r", 0.15);
    plt::arrow(end_x, end_y, cos(end_yaw), sin(end_yaw), "r", 0.15);
    plt::legend();
    plt::axis("equal");
    plt::grid(true);
    plt::show();
}

int main(int argc, char** argv)
{
    std::string mode = "static";
    if (argc > 1) {
        mode = argv[1];
    }

    if (mode == "static") {
        bezier_path();
    } else {
        dynamic_effect();
    }

    return 0;
}
