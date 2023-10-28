#include <cmath>
#include <random>
#include <tuple>
#include <vector>
#include <unistd.h>

#include <fmt/core.h>
#include <Eigen/Core>

#include "PathFinderController.hpp"
#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

PathFinderController controller = PathFinderController(9, 15, 3);
double dt = 0.01;

constexpr int MAX_LINEAR_SPEED = 15;
constexpr int MAX_ANGULAR_SPEED = 10;
bool show_animation = true;

void plot_vehicle(double x, double y, double theta, 
            const std::vector<double>& x_traj, const std::vector<double>& y_traj)
{
    Vector3d p1_i(0.5, 0, 1);
    Vector3d p2_i(-0.5, 0.25, 1);
    Vector3d p3_i(-0.5, -0.25, 1);

    Matrix3d T = Utils::transformation_matrix2d(x, y, theta);
    Vector3d p1 = T * p1_i;
    Vector3d p2 = T * p2_i;
    Vector3d p3 = T * p3_i;
    plt::plot({p1[0], p2[0]}, {p1[1], p2[1]}, "k-");
    plt::plot({p2[0], p3[0]}, {p2[1], p3[1]}, "k-");
    plt::plot({p1[0], p3[0]}, {p1[1], p3[1]}, "k-");
    plt::plot(x_traj, y_traj, "b--");
    plt::xlim(0, 20);
    plt::ylim(0, 20);
    plt::pause(dt);
}

void move_to_pose(double x_start, double y_start, double theta_start,
                double x_goal, double y_goal, double theta_goal)
{
    double x = x_start;
    double y = y_start;
    double theta = theta_start;
    double x_diff = x_goal - x;
    double y_diff = y_goal - y;

    std::vector<double> x_traj;
    std::vector<double> y_traj;

    double rho = hypot(x_diff, y_diff);
    int epoch = 0;
    while (rho > 0.005) {
        double v;
        double w;

        x_traj.emplace_back(x);
        y_traj.emplace_back(y);
        x_diff = x_goal - x;
        y_diff = y_goal - y;
        
        std::tie(rho, v, w) = controller.calc_control_command(x_diff, y_diff, theta, theta_goal);
        if (abs(v) > MAX_LINEAR_SPEED) {
            v = Utils::sign(v) * MAX_LINEAR_SPEED;
        }
        if (abs(w) > MAX_ANGULAR_SPEED) {
            w = Utils::sign(w) * MAX_ANGULAR_SPEED;
        }
        theta = theta + w * dt;
        x = x + v * cos(theta) * dt;
        y = y + v * sin(theta) * dt;

        if (show_animation) {
            plt::clf();
            plt::arrow(x_start, y_start, cos(theta_start), sin(theta_start), "r");
            plt::arrow(x_goal, y_goal, cos(theta_goal), sin(theta_goal), "g");
            plt::title("move_to_pose");
            plot_vehicle(x, y, theta, x_traj, y_traj);
        }
        epoch++;
        if (epoch > 320) {
            fmt::print("Planning failed, current deviation: {:.5f}\n", rho);
            break;
        }
    }
}

int main(int argv, char** argc)
{
    std::random_device seed;
    std::mt19937 eigine(seed());
    std::uniform_real_distribution<double> distrib(0, 1);

    for (int i = 0; i < 5; ++i) {
        double x_start = 18 * distrib(eigine) + 1.0;
        double y_start = 18 * distrib(eigine) + 1.0;
        double theta_start = 2 * M_PI * distrib(eigine) - M_PI;
        double x_goal = 18 * distrib(eigine) + 1.0;
        double y_goal = 18 * distrib(eigine) + 1.0;
        double theta_goal = 2 * M_PI * distrib(eigine) - M_PI;
        fmt::print("Initial x: {:.2f} m\tInitial y: {:.2f} m\tInitial theta: {:.2f} rad\n",
            x_start, y_start, theta_start);
        fmt::print("Goal x: {:.2f} m\t\tGoal y: {:.2f} m\t\tGoal theta: {:.2f} rad\n",
            x_goal, x_goal, theta_goal);
        move_to_pose(x_start, y_start, theta_start, x_goal, y_goal, theta_goal);
        fmt::print("\n");
    }
    

    return 0;
}

