#include <cmath>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <fmt/core.h>

#include "utils.hpp"
#include "matplotlibcpp.h"

using std::vector;
using std::string;
using namespace Eigen;
namespace plt = matplotlibcpp;

typedef Matrix<double, 5, 1> Vector5d;
constexpr bool show_animation = true;

double mod2pi(double theta)
{
    return std::fmod(theta, 2 * M_PI);
}

Vector5d calc_trig_funcs(double alpha, double beta)
{
    Vector5d trig;
    trig << sin(alpha), sin(beta), cos(alpha),
            cos(beta), cos(alpha - beta);
    
    return trig;
}

Vector3d lsl(double alpha, double beta, double d)
{
    Vector5d trig = calc_trig_funcs(alpha, beta);

    double p_squared = 2 + d * d - (2 * trig[4]) + (2 * d * (trig[0] - trig[1]));
    if (p_squared < 0) {
        return {-1, -1, -1};
    }
    double tmp = atan2((trig[3] - trig[2]), d + trig[0] - trig[1]);
    double d1 = mod2pi(-alpha + tmp);
    double d2 = sqrt(p_squared);
    double d3 = mod2pi(beta - tmp);
    
    return {d1, d2, d3};
}

Vector3d rsr(double alpha, double beta, double d)
{
    Vector5d trig = calc_trig_funcs(alpha, beta);

    double p_squared = 2 + d * d - (2 * trig[4]) + (2 * d * (trig[1] - trig[0]));
    if (p_squared < 0) {
        return {-1, -1, -1};
    }
    double tmp = atan2((trig[2] - trig[3]), d - trig[0] + trig[1]);
    double d1 = mod2pi(alpha - tmp);
    double d2 = sqrt(p_squared);
    double d3 = mod2pi(-beta + tmp);
    
    return {d1, d2, d3};
}

Vector3d lsr(double alpha, double beta, double d)
{
    Vector5d trig = calc_trig_funcs(alpha, beta);

    double p_squared = -2 + d * d + (2 * trig[4]) + (2 * d * (trig[0] + trig[1]));
    if (p_squared < 0) {
        return {-1, -1, -1};
    }
    double d1 = sqrt(p_squared);
    double tmp = atan2((-trig[2] - trig[3]), (d + trig[0] + trig[1])) - atan2(-2.0, d1);
    double d2 = mod2pi(-alpha + tmp);
    double d3 = mod2pi(-mod2pi(beta) + tmp);
    
    return {d2, d1, d3};
}

Vector3d rsl(double alpha, double beta, double d)
{
    Vector5d trig = calc_trig_funcs(alpha, beta);

    double p_squared = d * d - 2 + (2 * trig[4]) - (2 * d * (trig[0] + trig[1]));
    if (p_squared < 0) {
        return {-1, -1, -1};
    }
    double d1 = sqrt(p_squared);
    double tmp = atan2((trig[2] + trig[3]), (d - trig[0] - trig[1])) - atan2(2.0, d1);
    double d2 = mod2pi(alpha - tmp);
    double d3 = mod2pi(beta - tmp);

    return {d2, d1, d3};
}

Vector3d rlr(double alpha, double beta, double d)
{
    Vector5d trig = calc_trig_funcs(alpha, beta);

    double tmp = (6.0 - d * d + 2.0 * trig[4] + 2.0 * d * (trig[0] - trig[1])) / 8.0;
    if (tmp < 0) {
        return {-1, -1, -1};
    }
    double d2 = mod2pi(2 * M_PI - acos(tmp));
    double d1 = mod2pi(alpha - atan2(trig[2] - trig[3], d - trig[0] + trig[1]) + d2 / 2.0);
    double d3 = mod2pi(alpha - beta - d1 + d2);
    
    return {d1, d2, d3};
}

Vector3d lrl(double alpha, double beta, double d)
{
    Vector5d trig = calc_trig_funcs(alpha, beta);

    double tmp = (6.0 - d * d + 2.0 * trig[4] + 2.0 * d * (-trig[0] + trig[1])) / 8.0;
    if (tmp < 0) {
        return {-1, -1, -1};
    }
    double d2 = mod2pi(2 * M_PI - acos(tmp));
    double d1 = mod2pi(-alpha - atan2(trig[2] - trig[3], d + trig[0] - trig[1]) + d2 / 2.0);
    double d3 = mod2pi(mod2pi(beta) - alpha - d1 + mod2pi(d2));
    
    return {d1, d2, d3};
}

Vector3d interpolate(double length, char mode, double max_curvature, Vector3d origin)
{
    Vector3d inter;

    if (mode == 'S') {
        inter[0] = origin[0] + length / max_curvature * cos(origin[2]);
        inter[1] = origin[1] + length / max_curvature * sin(origin[2]);
        inter[2] = origin[2];
    } else {
        double ldx = sin(length) / max_curvature;
        double ldy = 0.0;
        if (mode == 'L') {
            ldy = (1.0 - cos(length)) / max_curvature;
        } else if (mode == 'R') {
            ldy = (1.0 - cos(length)) / -max_curvature;
        }
        double gdx = cos(-origin[2]) * ldx + sin(-origin[2]) * ldy;
        double gdy = -sin(-origin[2]) * ldx + cos(-origin[2]) * ldy;
        inter[0] = origin[0] + gdx;
        inter[1] = origin[1] + gdy;
        if (mode == 'L') {
            inter[2] = origin[2] + length;
        } else if (mode == 'R') {
            inter[2] = origin[2] - length;
        }
    }

    return inter;
}

vector<vector<double>> generate_local_course(
    Vector3d lengths, string modes, double max_curvature, double step_size)
{
    vector<vector<double>> course = {{0.}, {0.}, {0.}};
    
    for (size_t idx = 0; idx < 3; ++idx) {
        if (lengths[idx] == 0.0) {
            continue;
        }
        Vector3d origin = {course[0].back(), course[1].back(), course[2].back()};
        double current_length = step_size;
        while (abs(current_length + step_size) <= abs(lengths[idx])) {

            Vector3d inter = interpolate(
                    current_length, modes[idx], max_curvature, origin);
            course[0].push_back(inter[0]);
            course[1].push_back(inter[1]);
            course[2].push_back(inter[2]);
            current_length += step_size;
        }
        Vector3d inter = interpolate(
            lengths[idx], modes[idx], max_curvature, origin);
        course[0].push_back(inter[0]);
        course[1].push_back(inter[1]);
        course[2].push_back(inter[2]);
    }

    return course;
}

vector<vector<double>> dubins_path_planning_from_origin(
    double end_x, double end_y, double end_yaw, double curvature,
    double step_size, const vector<string>& planning_funcs, string& m)
{
    double dx = end_x;
    double dy = end_y;
    double d = hypot(dx, dy) * curvature;
    double theta = mod2pi(atan2(dy, dx));
    double alpha = mod2pi(-theta);
    double beta = mod2pi(end_yaw - theta);
    double best_cost = std::numeric_limits<double>::max();
    Vector3d best_dis;
    string best_mode;
    
    for (string planner : planning_funcs) {
        Vector3d distance;
        if (planner == "LSL") {
            distance = lsl(alpha, beta, d);
        } else if (planner == "RSR") {
            distance = rsr(alpha, beta, d);
        } else if (planner == "LSR") {
            distance = lsr(alpha, beta, d);
        } else if (planner == "RSL") {
            distance = rsl(alpha, beta, d);
        } else if (planner == "RLR") {
            distance = rlr(alpha, beta, d);
        } else if (planner == "LRL") {
            distance = lrl(alpha, beta, d);
        } else {
            fmt::print("Invalid mode !");
            continue;
        }
        if (distance[0] < 0 && distance[1] < 0 && distance[2] < 0) {
            continue;
        }

        double cost = (abs(distance[0]) + abs(distance[1]) + abs(distance[2]));
        if (cost < best_cost) {
            best_cost = cost;
            best_dis = distance;
            best_mode = planner;
        }
    }

    m = best_mode;
    vector<vector<double>> course = generate_local_course(
                    best_dis, best_mode, curvature, step_size);
    
    return course;
}

vector<vector<double>> convert_xy(const vector<vector<double>>& path,
                                        Matrix2d rot, Vector3d start)
{
    vector<vector<double>> convert_path(3);

    for (size_t idx = 0; idx < path[0].size(); ++idx) {
        Vector2d point(path[0][idx], path[1][idx]);
        Matrix<double, 1, 2> converted_xy = point.transpose() * rot;
        convert_path[0].emplace_back(converted_xy(0, 0) + start[0]);
        convert_path[1].emplace_back(converted_xy(0, 1) + start[1]);
        convert_path[2].emplace_back(path[2][idx] + start[2]);
    }

    return convert_path;
}

vector<vector<double>> plan_dubins_path(Vector3d start, Vector3d goal, 
        double curvature, string& mode, double step_size = 0.1,
        vector<string> selected_types = {})
{
    double s_x = start[0];
    double s_y = start[1];
    double s_yaw = start[2];
    double g_x = goal[0];
    double g_y = goal[1];
    double g_yaw = goal[2];

    if (selected_types.empty()) {
        selected_types = {"LSL", "RSR", "LSR", "RSL", "RLR", "LRL"};
    }

    Matrix2d l_rot = utils::rotation_matrix2d(s_yaw);
    Vector2d s_to_g(g_x - s_x, g_y - s_y);
    Matrix<double, 1, 2> le_xy = s_to_g.transpose() * l_rot;
    double local_goal_x = le_xy(0, 0);
    double local_goal_y = le_xy(0, 1);
    double local_goal_yaw = g_yaw - s_yaw;

    vector<vector<double>> path = dubins_path_planning_from_origin(
        local_goal_x, local_goal_y, local_goal_yaw, curvature,
        step_size, selected_types, mode);
    Matrix2d rot = utils::rotation_matrix2d(-s_yaw);
    
    return convert_xy(path, rot, start);
}

int main(int argc, char** argv)
{
    Vector3d start(1., 1., M_PI_4);
    Vector3d goal(-3., -3., -M_PI_4);
    double curvature = 1.0;
    string mode;
    utils::VehicleConfig vc(0.25);

    vector<vector<double>> path = plan_dubins_path(start, goal, curvature, mode);

    double steer = 0.0;
    if (show_animation) {
        for (size_t idx = 0; idx < path[0].size(); ++idx) {
            plt::cla();
            plt::named_plot(mode, path[0], path[1]);
            plt::arrow(start[0], start[1], cos(start[2]), sin(start[2]), "r", 0.075);
            plt::arrow(goal[0], goal[1], cos(goal[2]), sin(goal[2]), "g", 0.075);

            if (idx < path[0].size() - 2) {
                double dy = (path[2][idx + 1] - path[2][idx]) / 0.5;
                steer = -utils::pi_2_pi(atan(-3.5 * dy));
            } else {
                steer = 0.0;
            }
            utils::draw_vehicle({path[0][idx], path[1][idx], path[2][idx]}, steer, vc);
            plt::legend({{"loc", "upper right"}});
            plt::grid(true);
            plt::axis("equal");
            plt::title("Dubins Path Planning");
            plt::pause(0.001);
        }
        plt::show();
    }

    return 0;
}
