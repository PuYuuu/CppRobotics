#include <cmath>
#include <random>
#include <tuple>
#include <vector>
#include <string>
#include <unistd.h>

#include <fmt/core.h>
#include <Eigen/Core>

#include "PathFinderController.hpp"
#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using namespace Eigen;
using std::string;
namespace plt = matplotlibcpp;

constexpr int TIME_DURATION = 1000;
constexpr double TIME_STEP = 0.01;
constexpr double AT_TARGET_ACCEPTANCE_THRESHOLD = 0.01;
constexpr bool SHOW_ANIMATION = true;
constexpr int PLOT_WINDOW_SIZE_X = 20;
constexpr int PLOT_WINDOW_SIZE_Y = 20;
constexpr int PLOT_FONT_SIZE = 8;

bool simulation_running = true;
bool all_robots_are_at_target = false;

class Pose
{
public:
    double x;
    double y;
    double theta;
public:
    Pose() : x(0), y(0), theta(0) {}
    Pose(double _x, double _y, double _theta) :
        x(_x), y(_y), theta(_theta) {}
    ~Pose() {}
};

class Robot
{
private:
    double max_linear_speed;
    double max_angular_speed;
    PathFinderController path_finder_controller;
public:
    string name;
    string color;
    bool is_at_target;
    std::vector<double> x_traj;
    std::vector<double> y_traj;
    Pose pose;
    Pose pose_start;
    Pose pose_target;

    Robot(string _name, string _color, double _max_linear_speed, double _max_angular_speed,
                 PathFinderController _controller) {
        name = _name;
        color = _color;
        max_linear_speed = _max_linear_speed;
        max_angular_speed = _max_angular_speed;
        path_finder_controller = _controller;
        is_at_target = false;
    }
    ~Robot() {}

    void set_start_target_poses(Pose _pose_start, Pose _pose_target);
    void move(double dt);

};

void Robot::set_start_target_poses(Pose _pose_start, Pose _pose_target)
{
    pose_start = _pose_start;
    pose_target = _pose_target;
    pose = pose_start;
}

void Robot::move(double dt)
{
    double linear_velocity;
    double angular_velocity;
    double rho;

    x_traj.emplace_back(pose.x);
    y_traj.emplace_back(pose.y);
    std::tie(rho, linear_velocity, angular_velocity) = path_finder_controller.calc_control_command(
        pose_target.x - pose.x, pose_target.y - pose.y, pose.theta, pose_target.theta);
    if (rho < AT_TARGET_ACCEPTANCE_THRESHOLD) {
        is_at_target = true;
    }
    if (abs(linear_velocity) > max_linear_speed) {
        linear_velocity = Utils::sign(linear_velocity) * max_linear_speed;
    }
    if (abs(angular_velocity) > max_angular_speed) {
        angular_velocity = Utils::sign(angular_velocity) * max_angular_speed;
    }

    pose.theta = pose.theta + angular_velocity * dt;
    pose.x = pose.x + linear_velocity * cos(pose.theta) * dt;
    pose.y = pose.y + linear_velocity * sin(pose.theta) * dt;
}

void plot_vehicle(double x, double y, double theta, string color,
            const std::vector<double>& x_traj, const std::vector<double>& y_traj)
{
    Vector3d p1_i(0.5, 0, 1);
    Vector3d p2_i(-0.5, 0.25, 1);
    Vector3d p3_i(-0.5, -0.25, 1);

    Matrix3d T = Utils::transformation_matrix2d(x, y, theta);
    Vector3d p1 = T * p1_i;
    Vector3d p2 = T * p2_i;
    Vector3d p3 = T * p3_i;
    plt::plot({p1[0], p2[0]}, {p1[1], p2[1]}, color + "-");
    plt::plot({p2[0], p3[0]}, {p2[1], p3[1]}, color + "-");
    plt::plot({p1[0], p3[0]}, {p1[1], p3[1]}, color + "-");
    plt::plot(x_traj, y_traj, color + "--");
}

void run_simulation(std::vector<Robot>& robots)
{
    double time = 0;
    int at_target_robot_num = 0;
    std::vector<string> robots_names;

    for (const Robot& robo : robots) {
        robots_names.emplace_back(robo.name);
    }

    while (simulation_running && time < TIME_DURATION) {
        time += TIME_STEP;

        for (Robot& robo : robots) {
            if (!robo.is_at_target) {
                robo.move(TIME_STEP);
                if (robo.is_at_target) {
                    at_target_robot_num += 1;
                }
            }
        }
        if (at_target_robot_num == robots.size()) {
            simulation_running = false;
        }
        
        if (SHOW_ANIMATION) {
            plt::clf();
            plt::xlim(0, PLOT_WINDOW_SIZE_X);
            plt::ylim(0, PLOT_WINDOW_SIZE_Y);
            plt::text(0.3, double(PLOT_WINDOW_SIZE_Y - 1), fmt::format("Time: {:.2f}", time));
            plt::text(0.3, double(PLOT_WINDOW_SIZE_Y - 2), fmt::format("Reached target robot num: {}",
                                                                            at_target_robot_num));
            for (Robot& robo : robots) {
                plt::arrow(robo.pose_start.x, robo.pose_start.y, cos(robo.pose_start.theta), 
                    sin(robo.pose_start.theta), "r");
                plt::arrow(robo.pose_target.x, robo.pose_target.y, cos(robo.pose_target.theta), 
                    sin(robo.pose_target.theta), "g");
                plot_vehicle(robo.pose.x, robo.pose.y, robo.pose.theta, robo.color,
                    robo.x_traj, robo.y_traj);
            }
            plt::pause(TIME_STEP);
        }
    }
}

int main(int argc, char** argv)
{
    Pose pose_target = Pose(15, 15, -1);
    Pose pose_start_1 = Pose(5, 2, 0);
    Pose pose_start_2 = Pose(5, 2, 0);
    Pose pose_start_3 = Pose(5, 2, 0);

    PathFinderController controller_1 = PathFinderController(5, 8, 2);
    PathFinderController controller_2 = PathFinderController(5, 16, 4);
    PathFinderController controller_3 = PathFinderController(10, 25, 6);

    Robot robot_1 = Robot("Yellow Robot", "y", 12, 5, controller_1);
    Robot robot_2 = Robot("Black Robot", "k", 16, 5, controller_2);
    Robot robot_3 = Robot("Blue Robot", "b", 20, 5, controller_3);

    robot_1.set_start_target_poses(pose_start_1, pose_target);
    robot_2.set_start_target_poses(pose_start_2, pose_target);
    robot_3.set_start_target_poses(pose_start_3, pose_target);

    std::vector<Robot> robots = {robot_1, robot_2, robot_3};
    run_simulation(robots);

    return 0;
}
