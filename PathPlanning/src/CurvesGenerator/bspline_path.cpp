#include <fmt/core.h>

#include <Eigen/Core>
#include <cmath>
#include <string>
#include <vector>

#include "matplotlibcpp.h"
#include "utils.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

enum class Type { Uniform, QUniform };

// There are problems in calculating the first and second
// derivativesin non-equidistant sequences.
void plot_curvature(const vector<vector<double>>& path, std::string style = "-c",
                    std::string label = "Curvature") {
    size_t size = path[0].size();
    vector<double> first_order;
    // vector<double> second_order;
    vector<double> yaw;
    vector<double> curvature;

    for (size_t idx = 0; idx < size - 1; ++idx) {
        double dydx = (path[1][idx + 1] - path[1][idx]) / (path[0][idx + 1] - path[0][idx]);
        first_order.emplace_back(dydx);
        yaw.emplace_back(atan(dydx));
    }
    yaw.push_back(yaw.back());
    first_order.push_back(first_order.back());

    for (size_t idx = 1; idx < size - 1; ++idx) {
        double ddyddx = (path[1][idx + 1] - 2 * path[1][idx] + path[1][idx - 1]) /
                        ((path[0][idx + 1] - path[0][idx]) * (path[0][idx + 1] - path[0][idx]));
        // second_order.emplace_back(ddyddx);
        curvature.emplace_back(abs(ddyddx) / pow((1 + pow(first_order[idx], 2)), 1.5));
    }
    curvature.push_back(curvature.back());
    curvature.push_back(curvature.back());

    vector<vector<double>> cxy(2);
    for (size_t idx = 0; idx < size; ++idx) {
        cxy[0].push_back(path[0][idx] + 2 * curvature[idx] * cos(yaw[idx]));
        cxy[1].push_back(path[1][idx] + 2 * curvature[idx] * sin(yaw[idx]));
    }

    plt::named_plot(label, cxy[0], cxy[1], style);
    for (size_t idx = 0; idx < size; ++idx) {
        plt::plot({path[0][idx], cxy[0][idx]}, {path[1][idx], cxy[1][idx]}, style);
    }
}

double bspline_bfunc(int i, int k, double uu, const vector<double>& u) {
    double bfunc = 0.0;
    if (k == 1) {
        if (u[i] <= uu && uu < u[i + 1]) {
            bfunc = 1.0;
        } else {
            bfunc = 0.0;
        }
    } else if (k >= 2) {
        double A = 0.0;
        double B = 0.0;

        if (u[i + k - 1] - u[i] == 0.0) {
            A = 0.0;
        } else {
            A = (uu - u[i]) / (u[i + k - 1] - u[i]);
        }

        if (u[i + k] - u[i + 1] == 0.0) {
            B = 0.0;
        } else {
            B = (u[i + k] - uu) / (u[i + k] - u[i + 1]);
        }

        bfunc = A * bspline_bfunc(i, k - 1, uu, u) + B * bspline_bfunc(i + 1, k - 1, uu, u);
    }

    return bfunc;
}

vector<vector<double>> plan_Bspline_path(int k, Type _type, vector<vector<double>> p) {
    double delta_u = 0.01;
    int order = k;
    int n = p[0].size() - 1;
    if (order > n + 1 || p.empty()) {
        throw std::invalid_argument("Bspline illegal parameter !");
        return {};
    }

    vector<double> u;
    Type type = _type;
    vector<vector<double>> control_point = p;

    double u_tmp = 0.0;
    u.push_back(u_tmp);

    if (type == Type::Uniform) {
        double dis_u = 1.0 / (k + n);
        for (size_t i = 1; i < n + k + 1; ++i) {
            u_tmp += dis_u;
            u.push_back(u_tmp);
        }
    } else if (type == Type::QUniform) {
        int j = 3;
        double dis_u = 1.0 / (k + n - (j - 1) * 2);
        for (size_t i = 1; i < j; ++i) {
            u.push_back(u_tmp);
        }
        for (int i = j; i < n + k - j + 2; i++) {
            u_tmp += dis_u;
            u.push_back(u_tmp);
        }
        for (int i = n + k - j + 2; i < n + k + 1; i++) {
            u.push_back(u_tmp);
        }
    }

    double u_begin = u[k - 1];
    double u_end = u[n + 1];

    vector<vector<double>> path(2);
    for (double uu = u_begin; uu <= u_end; uu += delta_u) {
        Vector2d p_u(0.0, 0.0);
        for (size_t idx = 0; idx < n + 1; ++idx) {
            double xtmp = control_point[0][idx];
            double ytmp = control_point[1][idx];
            double bfunc_tmp = bspline_bfunc(idx, order, uu, u);
            p_u[0] += xtmp * bfunc_tmp;
            p_u[1] += ytmp * bfunc_tmp;
        }
        path[0].push_back(p_u[0]);
        path[1].push_back(p_u[1]);
    }

    return path;
}

int main(int argc, char** argv) {
    vector<vector<double>> way_point = {{0., 3., 6., 2., 1., 4.}, {0., -3., 0., 1., 3., 4.}};

    vector<vector<double>> path = plan_Bspline_path(3, Type::QUniform, way_point);

    plt::named_plot("Interpolated B-Spline path", path[0], path[1], "-b");
    plt::named_plot("way points", way_point[0], way_point[1], "-og");
    // plot_curvature(path);
    plt::title("B-Spline Interpolation");
    plt::legend();
    plt::xlabel("x[m]");
    plt::ylabel("y[m]");
    plt::grid(true);
    plt::axis("equal");
    plt::show();

    return 0;
}
