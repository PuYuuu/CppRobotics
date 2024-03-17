#pragma once

#ifndef __PathFinderController_HPP
#define __PathFinderController_HPP

#include <cmath>
#include <tuple>

class PathFinderController {
private:
    double kp_rho;
    double kp_alpha;
    double kp_beta;

public:
    PathFinderController() {}
    PathFinderController(double _kp_rho, double _kp_alpha, double _kp_beta) {
        kp_rho = _kp_rho;
        kp_alpha = _kp_alpha;
        kp_beta = _kp_beta;
    }
    ~PathFinderController() {}

    std::tuple<double, double, double> calc_control_command(double x_diff, double y_diff,
                                                            double theta, double theta_goal);
};

std::tuple<double, double, double> PathFinderController::calc_control_command(double x_diff,
                                                                              double y_diff,
                                                                              double theta,
                                                                              double theta_goal) {
    double rho = hypot(x_diff, y_diff);
    double alpha = std::fmod(std::atan2(y_diff, x_diff) - theta + M_PI, 2 * M_PI) - M_PI;
    double beta = std::fmod(theta_goal - theta - alpha + M_PI, 2 * M_PI) - M_PI;
    double v = kp_rho * rho;
    double w = kp_alpha * alpha - kp_beta * beta;
    if (alpha > M_PI_2 || alpha < -M_PI_2) {
        v = -v;
    }

    return std::make_tuple(rho, v, w);
}

#endif
