#include <fmt/core.h>

#include <Eigen/Eigen>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "PathPlanning/include/cubic_spline.hpp"
#include "matplotlibcpp.h"
#include "utils.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr double DT = 0.1;
constexpr double MAX_SIM_TIME = 100.0;
constexpr bool show_animation = true;

double k = 0.5;   // control gain
double Kp = 1.0;  // speed proportional gain

std::pair<size_t, double> calc_target_index(const utils::VehicleState& state,
                                            const vector<double>& cx, const vector<double>& cy) {
    double fx = state.x + (state.vc.RF) * cos(state.yaw);
    double fy = state.y + (state.vc.RF) * sin(state.yaw);
    double min_d = std::numeric_limits<double>::max();
    size_t target_idx = -1;
    Vector2d error_vec;

    for (size_t idx = 0; idx < cx.size(); ++idx) {
        double dx = fx - cx[idx];
        double dy = fy - cy[idx];
        double d = hypot(dx, dy);
        if (d <= min_d) {
            min_d = d;
            target_idx = idx;
            error_vec << dx, dy;
        }
    }

    Vector2d front_axle_vec(-cos(state.yaw + M_PI_2), -sin(state.yaw + M_PI_2));
    double error_front_axle = error_vec.dot(front_axle_vec);

    return std::make_pair(target_idx, error_front_axle);
}

double stanley_control(const utils::VehicleState& state, const vector<double>& cx,
                       const vector<double>& cy, const vector<double>& cyaw,
                       size_t& last_target_idx) {
    auto _target = calc_target_index(state, cx, cy);
    size_t current_target_idx = _target.first;
    double error_front_axle = _target.second;

    if (last_target_idx >= current_target_idx) {
        current_target_idx = last_target_idx;
    }

    // Sometimes you need to set a scaling factor for theta_e, e.g. 0.8
    double theta_e = utils::pi_2_pi(cyaw[current_target_idx] - state.yaw) * 0.8;
    double theta_d = atan2(k * error_front_axle, state.v);
    double delta = theta_e + theta_d;
    last_target_idx = current_target_idx;

    return delta;
}

int main(int argc, char** argv) {
    vector<double> ax = {0.0, 50.0, 50.0, 25.0, 30.0};
    vector<double> ay = {0.0, 0.0, -15.0, -10.0, 0.0};
    vector<vector<double>> c = CubicSpline2D::calc_spline_course(ax, ay, 0.1);
    vector<double> cx = c[0];
    vector<double> cy = c[1];
    vector<double> cyaw = c[2];
    double target_speed = 20.0 / 3.6;

    utils::VehicleConfig vc;
    utils::VehicleState state(vc, 0., 5., 20 * M_PI / 180.0, 0.);
    size_t last_idx = cx.size() - 1;
    double time = 0.0;
    vector<double> x = {state.x};
    vector<double> y = {state.y};
    vector<double> yaw = {state.yaw};
    vector<double> v = {state.v};
    vector<double> t = {0.0};
    auto _target = calc_target_index(state, cx, cy);
    size_t target_idx = _target.first;
    utils::TicToc t_m;

    while (MAX_SIM_TIME >= time && last_idx > target_idx) {
        double ai = Kp * (target_speed - state.v);
        double di = stanley_control(state, cx, cy, cyaw, target_idx);
        state.update(ai, di, DT);

        time += DT;
        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);
        t.push_back(time);

        if (show_animation) {
            plt::cla();
            plt::named_plot("course", cx, cy, "-r");
            plt::named_plot("trajectory", x, y, "-b");
            plt::plot({cx[target_idx]}, {cy[target_idx]}, "xg");
            utils::draw_vehicle({state.x, state.y, state.yaw}, di, vc);
            plt::legend();
            plt::axis("equal");
            plt::grid(true);
            plt::title("Stanley Controller Speed [km/h]:" +
                       std::to_string(state.v * 3.6).substr(0, 5));
            plt::pause(0.001);
        }
    }

    if (last_idx > target_idx) {
        fmt::print("Cannot reach goal\n", last_idx, target_idx);
        return 0;
    }

    fmt::print("stanley_controller run costtime: {:.3f} s\n", t_m.toc() / 1000);
    if (show_animation) {
        plt::cla();
        plt::named_plot("course", cx, cy, ".r");
        plt::named_plot("trajectory", x, y, "-b");
        plt::legend();
        plt::xlabel("x[m]");
        plt::ylabel("y[m]");
        plt::axis("equal");
        plt::grid(true);

        plt::figure();
        plt::plot(t, v, "-r");
        plt::xlabel("Time [s]");
        plt::ylabel("Speed [m/s]");
        plt::grid(true);
        plt::show();
    }

    return 0;
}
