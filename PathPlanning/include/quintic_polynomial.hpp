#pragma once
#ifndef __QUINTIC_POLYNOMIAL_HPP
#define __QUINTIC_POLYNOMIAL_HPP

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <cmath>

class QuinticPolynomial {
private:
    double a0;
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;

public:
    QuinticPolynomial(double xs, double vxs, double axs, double xe, double vxe, double axe,
                      double time) {
        a0 = xs;
        a1 = vxs;
        a2 = axs / 2.0;

        Eigen::Matrix3d A;
        A << pow(time, 3), pow(time, 4), pow(time, 5), 3 * pow(time, 2), 4 * pow(time, 3),
            5 * pow(time, 4), 6 * time, 12 * pow(time, 2), 20 * pow(time, 3);
        Eigen::Vector3d b;
        b << xe - a0 - a1 * time - a2 * pow(time, 2), vxe - a1 - 2 * a2 * time, axe - 2 * a2;

        Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
        a3 = x[0];
        a4 = x[1];
        a5 = x[2];
    }
    ~QuinticPolynomial() {}

    double calc_point(double t) {
        double xt = a0 + a1 * t + a2 * pow(t, 2) + a3 * pow(t, 3) + a4 * pow(t, 4) + a5 * pow(t, 5);

        return xt;
    }

    double calc_first_derivative(double t) {
        double xt = a1 + 2 * a2 * t + 3 * a3 * pow(t, 2) + 4 * a4 * pow(t, 3) + 5 * a5 * pow(t, 4);

        return xt;
    }

    double calc_second_derivative(double t) {
        double xt = 2 * a2 + 6 * a3 * t + 12 * a4 * pow(t, 2) + 20 * a5 * pow(t, 3);

        return xt;
    }

    double calc_third_derivative(double t) {
        double xt = 6 * a3 + 24 * a4 * t + 60 * a5 * pow(t, 2);

        return xt;
    }
};

#endif
