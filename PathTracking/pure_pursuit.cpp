#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <tuple>

#include <fmt/core.h>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::vector;
using std::tuple;
namespace plt = matplotlibcpp;
constexpr double WB = 2.9;
constexpr double dt = 0.1;
constexpr double k = 0.1;
constexpr double Lfc = 2.0;
constexpr double Kp = 1.0;
constexpr bool show_animation = true;

class RobotState
{
public:
    double x;
    double y;
    double yaw;
    double velocity;
    double rear_x;
    double rear_y;

    RobotState(double _x, double _y, double _yaw, double _v) : x(_x), y(_y), yaw(_yaw), velocity(_v) {
        rear_x = x - ((WB / 2) * cos(yaw));
        rear_y = y - ((WB / 2) * sin(yaw));
    }
    ~RobotState() {}

    void update(double acc, double delta);
    double calc_distance(double point_x, double point_y);
};

void RobotState::update(double acc, double delta)
{
    x += velocity * cos(yaw) * dt;
    y += velocity * sin(yaw) * dt;
    yaw += velocity / WB * tan(delta) * dt;
    velocity += acc * dt;
    rear_x = x - ((WB / 2) * cos(yaw));
    rear_y = y - ((WB / 2) * sin(yaw));
}

double RobotState::calc_distance(double point_x, double point_y)
{
    double dx = rear_x - point_x;
    double dy = rear_y - point_y;
    
    return hypot(dx, dy);
}     

class RobotStates
{
public:
    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    vector<double> velocity;
    vector<double> t;

    RobotStates() {}
    ~RobotStates() {}
    void append(double _t, RobotState state) {
        x.emplace_back(state.x);
        y.emplace_back(state.y);
        yaw.emplace_back(state.yaw);
        velocity.emplace_back(state.velocity);
        t.emplace_back(_t);
    }
};

class TargetCourse
{
public:
    vector<double> cx;
    vector<double> cy;
    int old_nearest_point_index;

    TargetCourse(vector<double> _cx, vector<double> _cy) : cx(_cx), cy(_cy) {
        old_nearest_point_index = -1;
    }
    ~TargetCourse() {}
    tuple<int, double> search_target_index(RobotState state);
};

tuple<int, double> TargetCourse::search_target_index(RobotState state)
{
    int ind = -1;

    if (old_nearest_point_index == -1) {
        double min_distance = std::numeric_limits<double>::max();
        for (int idx = 0; idx < cx.size(); ++idx) {
            double dx = state.rear_x - cx[idx];
            double dy = state.rear_y - cx[idx];
            double dis = hypot(dx, dy);
            if (dis < min_distance) {
                min_distance = dis;
                old_nearest_point_index = idx;
            }
        }
        ind = old_nearest_point_index;
    } else {
        ind = old_nearest_point_index;
        double distance_this_index = state.calc_distance(cx[ind], cy[ind]);
        while (true) {
            double distance_next_index = state.calc_distance(cx[ind + 1], cy[ind + 1]);
            if (distance_this_index < distance_next_index) {
                break;
            }
            ind = ((ind + 1) < cx.size()) ? (ind + 1) : ind;
            distance_this_index = distance_next_index;
        }
        old_nearest_point_index = ind;
    }

    double Lf = k * state.velocity + Lfc;
    while (Lf > state.calc_distance(cx[ind], cy[ind])) {
        if ((ind + 1) >= cx.size()) {
            break;
        }
        ind += 1;
    }

    return std::make_tuple(ind, Lf);
}

double proportional_control(double target, double current)
{
    return Kp * (target - current);
}

tuple<int, double> pure_pursuit_steer_control(RobotState& state, TargetCourse& trajectory, int pind)
{
    tuple<int, double> ret = trajectory.search_target_index(state);
    int ind = std::get<0>(ret);
    double Lf = std::get<1>(ret);

    if (pind >= ind) {
        ind = pind;
    }

    double tx;
    double ty;
    if (ind < trajectory.cx.size()) {
        tx = trajectory.cx[ind];
        ty = trajectory.cy[ind];
    } else {
        tx = trajectory.cx.back();
        ty = trajectory.cy.back();
        ind = trajectory.cx.size() - 1;
    }

    double alpha = atan2(ty - state.rear_y, tx - state.rear_x) - state.yaw;
    double delta = atan2(2.0 * WB * sin(alpha) / Lf, 1.0);

    return std::make_tuple(ind, delta);
}

int main(int argc, char** argv)
{
    Utils::TicToc t_m;
    vector<double> cx;
    vector<double> cy;
    for (double idx = 0; idx <= 50; idx += 0.5) {
        double tmp = sin(idx / 5.0) * idx;
        cx.emplace_back(idx);
        cy.emplace_back(tmp);
    }

    double target_speed = 10.0 / 3.6;
    double T = 100.0;
    RobotState state = RobotState(-0.0, -3.0, 0.0, 0.0);
    int lastIndex = cx.size() - 1;
    double time = 0.0;
    RobotStates states = RobotStates();
    states.append(time, state);
    TargetCourse target_course = TargetCourse(cx, cy);
    tuple<int, double> result = target_course.search_target_index(state);
    int target_id = std::get<0>(result);

    while (T >= time && lastIndex > target_id) {
        double ai = proportional_control(target_speed, state.velocity);
        result = pure_pursuit_steer_control(state, target_course, target_id);
        target_id = std::get<0>(result);
        double di = std::get<1>(result);

        state.update(ai, di);
        time += dt;
        states.append(time, state);

        if (show_animation) {
            plt::cla();
            plt::arrow(state.x, state.y, cos( state.yaw), sin(state.yaw), "k", 0.25);
            plt::named_plot("course", cx, cy, "-r");
            plt::named_plot("trajectory", states.x, states.y, "-b");
            plt::plot({cx[target_id]}, {cy[target_id]}, "xg");
            plt::axis("equal");
            plt::grid(true);
            plt::legend();
            // plt::xlim(-10, 55);
            plt::title("Speed[km/h]:" + std::to_string(state.velocity * 3.6).substr(0,4));
            plt::pause(0.001);
        }
    }

    fmt::print("pure_pursuit run costtime: {:.3f} s\n", t_m.toc() / 1000);
    if (show_animation) {
        plt::cla();
        plt::named_plot("course", cx, cy, ".r");
        plt::named_plot("trajectory", states.x, states.y, "-b");
        plt::legend();
        plt::xlabel("x[m]");
        plt::ylabel("y[m]");
        plt::axis("equal");
        plt::grid(true);
        plt::show();
    }

    return 0;
}
