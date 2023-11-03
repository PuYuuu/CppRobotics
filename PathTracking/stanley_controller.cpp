#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include <fmt/core.h>
#include <Eigen/Eigen>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"
#include "PathPlanning/CurvesGenerator/cubic_spline.hpp"

using std::vector;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;
constexpr double max_simulation_time = 100.0;

double k = 0.5;     // control gain
double Kp = 1.0;    // speed proportional gain
double dt = 0.1;    // [s] time difference
double L = 2.5;     // [m] Wheel base of vehicle
double max_steer = M_PI_2 / 3.0;  // [rad] max steering angle

vector<vector<double>> calc_spline_course(vector<double> x, vector<double> y, double ds = 0.1)
{
    CubicSpline2D sp(x, y);
    // rx, ry, ryaw
    vector<vector<double>> output(3);

    for (double s = sp.s.front(); s < sp.s.back(); s += ds) {
        Vector2d ixy = sp.calc_position(s);
        output[0].push_back(ixy[0]);
        output[1].push_back(ixy[1]);
        output[2].push_back(sp.calc_yaw(s));
    }

    return output;
}

double normalize_angle(double angle)
{
    while (angle > M_PI) {
        angle -= (2.0 * M_PI);
    }

    while (angle < -M_PI) {
        angle += (2.0 * M_PI);
    }

    return angle;
}

class RobotState
{ 
public:
    double x;
    double y;
    double yaw;
    double v;

    RobotState(double _x = 0, double _y = 0, double _yaw = 0, double _v = 0) :
        x(_x), y(_y), yaw(_yaw), v(_v) { }
    ~RobotState() {}
    void update(double acc, double delta);
};

void RobotState::update(double acc, double delta)
{
    if (delta < -max_steer) {
        delta = -max_steer;
    } else if (delta > max_steer) {
        delta = max_steer;
    }
    x += (v * cos(yaw) * dt);
    y += (v * sin(yaw) * dt);
    yaw += (v / L * tan(delta) * dt);
    yaw = normalize_angle(yaw);
    v += (acc * dt);
}

std::pair<size_t, double> calc_target_index(RobotState& state, vector<double> cx, vector<double> cy)
{
    double fx = state.x + L * cos(state.yaw);
    double fy = state.y + L * sin(state.yaw);
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

double stanley_control(RobotState& state, vector<double> cx, vector<double> cy,
    vector<double> cyaw, size_t& last_target_idx)
{
    auto _target = calc_target_index(state, cx, cy);
    size_t current_target_idx = _target.first;
    double error_front_axle = _target.second;

    if (last_target_idx >= current_target_idx) {
        current_target_idx = last_target_idx;
    }

    double theta_e = normalize_angle(cyaw[current_target_idx] - state.yaw);
    double theta_d = atan2(k * error_front_axle, state.v);
    double delta = theta_e + theta_d;
    last_target_idx = current_target_idx;

    return delta;
}

int main(int argc, char** argv)
{
    vector<double> ax = {0.0, 50.0, 50.0, 25.0, 30.0};
    vector<double> ay = {0.0, 0.0, -15.0, -10.0, 0.0};
    vector<vector<double>> c = calc_spline_course(ax, ay, 0.1);
    vector<double> cx = c[0];
    vector<double> cy = c[1];
    vector<double> cyaw = c[2];
    double target_speed = 20.0 / 3.6;
    
    RobotState state(-0.0, 5.0, 20 * M_PI / 180.0, 0.0);
    size_t last_idx = cx.size() - 1;
    double time = 0.0;
    vector<double> x = {state.x};
    vector<double> y = {state.y};
    vector<double> yaw = {state.yaw};
    vector<double> v = {state.v};
    vector<double> t = {0.0};
    auto _target = calc_target_index(state, cx, cy);
    size_t target_idx = _target.first;
    utils::VehicleConfig vc;
    utils::TicToc t_m;

    while (max_simulation_time >= time && last_idx > target_idx) {
        double ai = Kp * (target_speed - state.v);
        double di = stanley_control(state, cx, cy, cyaw, target_idx);
        state.update(ai, di);

        time += dt;
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
            plt::title("Stanley Controller Speed[km/h]:" + std::to_string(state.v * 3.6).substr(0, 5));
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
        plt::xlabel("Time[s]");
        plt::ylabel("Speed[m/s]");
        plt::grid(true);
        plt::show();
    }

    return 0;
}
