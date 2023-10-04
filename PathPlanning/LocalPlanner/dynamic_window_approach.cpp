#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>

#include <fmt/core.h>
#include <Eigen/Core>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::vector;
using std::shared_ptr;
using namespace Eigen;
// x(m), y(m), yaw(rad), v(m/s), omega(rad/s)
typedef Eigen::Matrix<double, 5, 1> RobotState;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

enum RobotType {Circle, Rectangle};

class Config
{
private:
    static Config instance;
public:
    double max_speed = 1.0;
    double min_speed = -0.5;
    double max_yaw_rate = 40.0 * M_PI / 180.0;
    double max_accel = 0.2;
    double max_delta_yaw_rate = 40.0 * M_PI / 180.0;
    double v_resolution = 0.05;
    double yaw_rate_resolution = 0.5 * M_PI / 180.0;
    double dt = 0.1;
    double predict_time = 3.0;
    double to_goal_cost_gain = 0.15;
    double speed_cost_gain = 1.0;
    double obstacle_cost_gain = 1.0;
    double robot_stuck_flag_cons = 0.001;
    RobotType robot_type = Rectangle;

    double robot_radius = 1.0;
    double robot_width = 0.5;
    double robot_length = 1.2;

    Config() {}
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
    ~Config() {}

    static Config* getInstance() {		
		return &instance;
	}
};
Config Config::instance;

Vector4d calc_dynamic_window(RobotState x, Config* config)
{
    Vector4d Vs(config->min_speed, config->max_speed,
                -config->max_yaw_rate, config->max_yaw_rate);

    Vector4d Vd(x[3] - config->max_accel * config->dt,
                x[3] + config->max_accel * config->dt,
                x[4] - config->max_delta_yaw_rate * config->dt,
                x[4] + config->max_delta_yaw_rate * config->dt);

    Vector4d dw(std::max(Vs[0], Vd[0]), std::min(Vs[1], Vd[1]),
                std::max(Vs[2], Vd[2]), std::min(Vs[3], Vd[3]));

    return Vs;
}

RobotState motion(RobotState x, double v, double w, double dt)
{
    RobotState x_pre;
    x_pre[2] = x[2] + w * dt;
    x_pre[0] = x[0] + v * cos(x[2]) * dt;
    x_pre[1] = x[1] + v * sin(x[2]) * dt;
    x_pre[3] = v;
    x_pre[4] = w;

    return x_pre;
}

vector<RobotState> predict_trajectory(RobotState x_init, double v, double w, Config* config)
{
    RobotState x = x_init;
    vector<RobotState> trajectory = {x_init};
    double time = 0.0;

    while (time <= config->predict_time) {
        x = motion(x, v, w, config->dt);
        trajectory.emplace_back(x);
        time += config->dt;
    }

    return trajectory;
}

double calc_to_goal_cost(vector<RobotState> trajectory, Vector2d goal)
{
    double dx = goal[0] - trajectory.back()[0];
    double dy = goal[1] - trajectory.back()[1];
    double error_angle = atan2(dy, dx);
    double cost_angle = error_angle - trajectory.back()[2];
    double cost = abs(atan2(sin(cost_angle), cos(cost_angle)));

    return cost;
}

double calc_obstacle_cost(vector<RobotState> trajectory, const vector<vector<double>>& ob, Config* config)
{
    double minr = std::numeric_limits<double>::max();

    for (int i = 0; i < ob.size(); ++i) {
        for (int j = 0; j < trajectory.size(); ++j) {
            double ox = ob[i][0];
            double oy = ob[i][1];
            double yaw = trajectory[j][2];
            double dx = trajectory[j][0] - ox;
            double dy = trajectory[j][1] - oy;
            double r = hypot(dx, dy);

            if (config->robot_type == Rectangle) {
                double obsx = ox * cos(yaw) - oy * sin(yaw);
                double obsy = ox * sin(yaw) + oy * cos(yaw);
                bool upper_check = (obsy <= config->robot_length / 2);
                bool right_check = (obsx <= config->robot_width / 2);
                bool bottom_check = (obsy >= -config->robot_length / 2);
                bool left_check = (obsx >= -config->robot_width / 2);
                if (upper_check && right_check && bottom_check && left_check) {
                    return std::numeric_limits<double>::max() / 2.0;
                }
            } else if (config->robot_type == Circle) {
                if (r <= config->robot_radius) {
                    return std::numeric_limits<double>::max() / 2.0;
                }
            }
            minr = std::min(minr, r);
        }
    }

    return 1.0 / minr;
}

vector<RobotState> calc_control_and_trajectory(RobotState x, Vector2d& control, Vector4d dw, 
                                Config* config, Vector2d goal, const vector<vector<double>>& ob)
{
    RobotState x_init = x;
    double min_cost = std::numeric_limits<double>::max();
    Vector2d best_u(0.0, 0.0);
    vector<RobotState> best_trajectory = {x};

    for (double v = dw[0]; v <= dw[1]; v += config->v_resolution) {
        for (double w = dw[2]; w <= dw[3]; w += config->yaw_rate_resolution) {
            
            vector<RobotState> trajectory = predict_trajectory(x_init, v, w, config);

            double to_goal_cost = config->to_goal_cost_gain * calc_to_goal_cost(trajectory, goal);
            double speed_cost = config->speed_cost_gain * (config->max_speed - trajectory.back()[3]);
            double ob_cost = config->obstacle_cost_gain * calc_obstacle_cost(trajectory, ob, config);
            double final_cost = to_goal_cost + speed_cost + ob_cost; 

            if (min_cost >= final_cost) {
                min_cost = final_cost;
                best_u << v, w;
                best_trajectory = trajectory;
                if (best_u[0] < config->robot_stuck_flag_cons && abs(x[3]) < config->robot_stuck_flag_cons) {
                    best_u[1] = -config->max_delta_yaw_rate;
                }
            }
        }
    }
    control = best_u;

    return best_trajectory;
}

vector<RobotState> dwa_control(RobotState x, Vector2d& control, Config* config, 
                            Vector2d goal, const vector<vector<double>>& ob)
{
    Vector4d dw = calc_dynamic_window(x, config);
    vector<RobotState> traj = calc_control_and_trajectory(x, control, dw, config, goal, ob);

    return traj;
}

void plot_robot(double x, double y, double yaw, Config* config)
{
    if (config->robot_type == Rectangle) {
        vector<Vector3d> outline(4);
        outline[0] << -config->robot_length / 2, config->robot_width / 2, 1;
        outline[1] << config->robot_length / 2, config->robot_width / 2, 1;
        outline[2] << config->robot_length / 2, -config->robot_width / 2, 1;
        outline[3] << -config->robot_length / 2, -config->robot_width / 2, 1;

        Matrix3d T = Utils::transformation_matrix2d(x, y, yaw);
        Vector3d p1 = T * outline[0];
        Vector3d p2 = T * outline[1];
        Vector3d p3 = T * outline[2];
        Vector3d p4 = T * outline[3];
        plt::plot({p1[0], p2[0]}, {p1[1], p2[1]}, "k-");
        plt::plot({p2[0], p3[0]}, {p2[1], p3[1]}, "k-");
        plt::plot({p3[0], p4[0]}, {p3[1], p4[1]}, "k-");
        plt::plot({p4[0], p1[0]}, {p4[1], p1[1]}, "k-");
    } else if (config->robot_type == Circle) {
        // todo
        // I'm confused how to draw a circle of specific size in matplotlibcpp
    }
}

int main (int argc, char** argv)
{
    RobotState x;
    x << 0., 0., 0., M_PI_4 / 2, 0., 0.;
    Vector2d goal(10, 10);
    Config* config = Config::getInstance();
    vector<vector<double>> obs = {{-1, -1}, {0, 2}, {4.0, 2.0}, {5.0, 4.0},
                            {5.0, 5.0}, {5.0, 6.0}, {5.0, 8.0}, {5.0, 9.0}, 
                            {8.0, 9.0}, {7.0, 9.0}, {8.0, 10.0}, {9.0, 11.0}, 
                            {12.0, 13.0}, {12.0, 12.0}, {15.0, 15.0}, {13.0, 13.0}};
    vector<RobotState> trajectory = {x};

    while (true) {
        Vector2d u;
        vector<RobotState> predicted_trajectory = dwa_control(x, u, config, goal, obs);
        x = motion(x, u[0], u[1], config->dt);

        if (show_animation) {
            plt::cla();
            vector<double> trajx, trajy;
            for (RobotState traj : predicted_trajectory) {
                trajx.emplace_back(traj[0]);
                trajy.emplace_back(traj[1]);
            }
            plt::plot(trajx, trajy, "-g");

            plt::plot({x[0]}, {x[1]}, "xr");
            plt::plot({goal[0]}, {goal[1]}, "xb");
            for (vector<double> ob : obs) {
                plt::plot({ob[0]}, {ob[1]}, "ok");
            }
            plot_robot(x[0], x[1], x[2], config);
            plt::axis("equal");
            plt::grid(true);
            plt::pause(0.0001);
        }
        double dist_to_goal = hypot(x[0] - goal[0], x[1] - goal[1]);
        if (dist_to_goal <= config->robot_radius) {
            break;
        }

        trajectory.push_back(x);
    }

    if (show_animation) {
        vector<double> trajx, trajy;
        for (RobotState traj : trajectory) {
            trajx.emplace_back(traj[0]);
            trajy.emplace_back(traj[1]);
        }
        plt::plot(trajx, trajy, "-r");
        plt::show();
    }

    return 0;
}
