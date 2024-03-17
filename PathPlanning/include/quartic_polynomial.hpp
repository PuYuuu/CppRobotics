#pragma once
#ifndef __QUARTIC_POLYNOMIAL_HPP
#define __QUARTIC_POLYNOMIAL_HPP

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <cmath>

class QuarticPolynomial {
private:
    double a0;
    double a1;
    double a2;
    double a3;
    double a4;

public:
    QuarticPolynomial(double xs, double vxs, double axs, double vxe, double axe, double time) {
        a0 = xs;
        a1 = vxs;
        a2 = axs / 2.0;

        Eigen::Matrix2d A;
        A << 3 * pow(time, 2), 4 * pow(time, 3), 6 * time, 12 * pow(time, 2);
        Eigen::Vector2d b;
        b << vxe - a1 - 2 * a2 * time, axe - 2 * a2;

        Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
        a3 = x[0];
        a4 = x[1];
    }
    ~QuarticPolynomial() {}

    double calc_point(double t) {
        double xt = a0 + a1 * t + a2 * pow(t, 2) + a3 * pow(t, 3) + a4 * pow(t, 4);

        return xt;
    }

    double calc_first_derivative(double t) {
        double xt = a1 + 2 * a2 * t + 3 * a3 * pow(t, 2) + 4 * a4 * pow(t, 3);

        return xt;
    }

    double calc_second_derivative(double t) {
        double xt = 2 * a2 + 6 * a3 * t + 12 * a4 * pow(t, 2);

        return xt;
    }

    double calc_third_derivative(double t) {
        double xt = 6 * a3 + 24 * a4 * t;

        return xt;
    }
};

#endif
