#include <cmath>
#include <vector>
#include <algorithm>

#include "utils.hpp"
#include "matplotlibcpp.h"
#include "cubic_spline.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

int main(int argc, char** argv)
{
    int mode = 0;
    if (argc > 1) {
        mode = atoi(argv[1]);
    }

    if (mode == 0) {
        vector<double> x = {0, 1, 2, 3, 4};
        vector<double> y = {1.7, -6, 5, 6.5, 0.0};
        CubicSpline sp(x, y);
        vector<double> xi;
        vector<double> yi;

        for (double i = 0; i < 4.001; i += 0.0999999999) {
            xi.push_back(i);
            yi.push_back(sp.calc_position(i));
        }
        plt::named_plot("Data points", x, y, "xb");
        plt::named_plot("Cubic spline interpolation", xi, yi, "r");
        plt::grid(true);
        plt::legend();
        plt::title("Cubic Spline");
        plt::show();
    } else {
        vector<double> x = {-2.5, 0.0, 2.5, 5.0, 7.5, 3.0, -1.0};
        vector<double> y = {0.7, -6, 5, 6.5, 0.0, 5.0, -2.0};
        CubicSpline2D sp(x, y);
        vector<double> rx;
        vector<double> ry;
        vector<double> ryaw;
        vector<double> rk;
        vector<double> s;

        for (double ds = 0; ds < sp.s.back(); ds += 0.1) {
            Vector2d ixy = sp.calc_position(ds);
            rx.push_back(ixy[0]);
            ry.push_back(ixy[1]);
            ryaw.push_back(sp.calc_yaw(ds));
            rk.push_back(sp.calc_curvature(ds));
            s.push_back(ds);
        }

        plt::figure();
        plt::named_plot("Data points", x, y, "xb");
        plt::named_plot("Cubic spline path", rx, ry, "-r");
        plt::grid(true);
        plt::axis("equal");
        plt::xlabel("x[m]");
        plt::ylabel("y[m]");
        plt::legend();
        plt::title("Cubic Spline Interpolation");

        plt::figure();
        plt::named_plot("yaw", s, ryaw, "-r");
        plt::grid(true);
        plt::legend();
        plt::xlabel("line length[m]");
        plt::ylabel("yaw angle[rad]");

        plt::figure();
        plt::named_plot("curvature", s, rk, "-r");
        plt::grid(true);
        plt::legend();
        plt::xlabel("line length[m]");
        plt::ylabel("curvature [1/m]");

        plt::show();
    }

    return 0;
}
