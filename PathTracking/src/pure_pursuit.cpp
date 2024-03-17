#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

#include "matplotlibcpp.h"
#include "utils.hpp"

using std::tuple;
using std::vector;
namespace plt = matplotlibcpp;
constexpr double DT = 0.1;
constexpr double k = 0.1;
constexpr double Lfc = 2.0;
constexpr double Kp = 1.0;
constexpr bool show_animation = true;

class TargetCourse {
public:
    vector<double> cx;
    vector<double> cy;
    int old_nearest_point_index;

    TargetCourse(vector<double> _cx, vector<double> _cy) : cx(_cx), cy(_cy) {
        old_nearest_point_index = -1;
    }
    ~TargetCourse() {}
    tuple<int, double> search_target_index(utils::VehicleState state);
};

tuple<int, double> TargetCourse::search_target_index(utils::VehicleState state) {
    int ind = -1;

    if (old_nearest_point_index == -1) {
        double min_distance = std::numeric_limits<double>::max();
        for (int idx = 0; idx < cx.size(); ++idx) {
            double dx = state.x - cx[idx];
            double dy = state.y - cx[idx];
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

    double Lf = k * state.v + Lfc;
    while (Lf > state.calc_distance(cx[ind], cy[ind])) {
        if ((ind + 1) >= cx.size()) {
            break;
        }
        ind += 1;
    }

    return std::make_tuple(ind, Lf);
}

double proportional_control(double target, double current) { return Kp * (target - current); }

tuple<int, double> pure_pursuit_steer_control(utils::VehicleState& state, TargetCourse& trajectory,
                                              int pind) {
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

    double alpha = atan2(ty - state.y, tx - state.x) - state.yaw;
    double delta = atan2(2.0 * state.vc.WB * sin(alpha) / Lf, 1.0);

    return std::make_tuple(ind, delta);
}

int main(int argc, char** argv) {
    vector<double> cx;
    vector<double> cy;
    for (double idx = 0; idx <= 50; idx += 0.5) {
        double tmp = sin(idx / 5.0) * idx / 3;
        cx.emplace_back(idx);
        cy.emplace_back(tmp);
    }

    double target_speed = 20.0 / 3.6;
    double T = 100.0;
    utils::VehicleConfig vc(0.8);
    utils::VehicleState state(vc, -0.0, -3.0, 0.0, 0.0);
    int lastIndex = cx.size() - 1;
    double time = 0.0;

    TargetCourse target_course = TargetCourse(cx, cy);
    tuple<int, double> result = target_course.search_target_index(state);
    int target_id = std::get<0>(result);
    utils::TicToc t_m;
    vector<double> x = {state.x};
    vector<double> y = {state.y};
    vector<double> yaw = {state.yaw};

    while (T >= time && lastIndex > target_id) {
        double ai = proportional_control(target_speed, state.v);
        result = pure_pursuit_steer_control(state, target_course, target_id);
        target_id = std::get<0>(result);
        double di = std::get<1>(result);

        state.update(ai, di, DT);
        time += DT;

        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);

        if (show_animation) {
            plt::cla();
            plt::named_plot("course", cx, cy, "-r");
            plt::named_plot("trajectory", x, y, "-b");
            plt::plot({cx[target_id]}, {cy[target_id]}, "xg");
            utils::draw_vehicle({state.x, state.y, state.yaw}, di, vc);
            plt::axis("equal");
            plt::grid(true);
            plt::legend();
            plt::title("Pure Pursuit Speed[km/h]:" + std::to_string(state.v * 3.6).substr(0, 4));
            plt::pause(0.001);
        }
    }

    fmt::print("pure_pursuit run costtime: {:.3f} s\n", t_m.toc() / 1000);
    if (show_animation) {
        plt::cla();
        plt::named_plot("course", cx, cy, ".r");
        plt::named_plot("trajectory", x, y, "-b");
        plt::legend();
        plt::xlabel("x[m]");
        plt::ylabel("y[m]");
        plt::axis("equal");
        plt::grid(true);
        plt::show();
    }

    return 0;
}
