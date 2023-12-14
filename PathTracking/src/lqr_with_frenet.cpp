#include <cmath>
#include <limits>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <fmt/core.h>

#include "utils.hpp"
#include "matplotlibcpp.h"
#include "PathPlanning/include/cubic_spline.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr double MAX_SIM_TIME = 500.0;
constexpr double GOAL_DIS = 0.3;
constexpr double DT = 0.1;
constexpr bool show_animation = true;

class TrajectoryAnalyzer
{
private:
    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    vector<double> k;
    size_t ind_old;
    size_t ind_end;
public:
    TrajectoryAnalyzer() {}
    TrajectoryAnalyzer(vector<double> _x, vector<double> _y, vector<double> _yaw, vector<double> _k)
        : x(_x), y(_y), yaw(_yaw), k(_k), ind_old(0), ind_end(_x.size() - 1) {}
    ~TrajectoryAnalyzer() {}

    Vector4d to_trajectory_frame(const utils::VehicleState& state);
};

Vector4d TrajectoryAnalyzer::to_trajectory_frame(const utils::VehicleState& vehicle_state)
{
    double x_cg = vehicle_state.x;
    double y_cg = vehicle_state.y;
    double cyaw = vehicle_state.yaw;
    // theta_e, e_cg, yaw_ref, k_ref
    Vector4d ret(0, 0, 0, 0);

    Vector3d min_dist(0, 0, std::numeric_limits<double>::max());
    size_t min_idx = ind_old;
    for (size_t idx = ind_old + 1; idx <= ind_end; ++idx) {
        double dx = x_cg - x[idx];
        double dy = y_cg - y[idx];
        double dist = hypot(dx, dy);
        if (dist < min_dist[2]) {
            min_dist << dx, dy, dist;
            min_idx = idx;
        }
    }
    ind_old = min_idx;

    Vector2d vec_axle_rot_90(cos(cyaw + M_PI_2), sin(cyaw + M_PI_2));
    Vector2d vec_path_2_cg(min_dist[0], min_dist[1]);
    if (vec_axle_rot_90.transpose() * vec_path_2_cg > 0) {
        ret[1] = min_dist[2];
    } else {
        ret[1] = -1 * min_dist[2];
    }

    ret[2] = yaw[ind_old];
    ret[0] = utils::pi_2_pi(cyaw - ret[2]);
    ret[3] = k[ind_old];

    return ret;
}

class LatController
{
private:
    double dt;
    double e_cg_old;
    double theta_e_old;

    MatrixXd A;
    MatrixXd B;
    MatrixXd Q;
    MatrixXd R;
public:
    explicit LatController(double _dt = 0.1) : dt(_dt), e_cg_old(0), theta_e_old(0) {
        A = Matrix4d::Zero();
        A(0, 0) = 1.0;
        A(0, 1) = dt;
        A(2, 2) = 1.0;
        A(2, 3) = dt;
        B = Matrix<double, 4, 1>::Zero();
        Q = Matrix4d::Identity();
        Q(0, 0) = 0.1;
        Q(2, 2) = 0.1;
        R = Matrix<double, 1, 1>::Identity();
    }
    LatController() = delete;
    ~LatController() {}

    double compute_input(
        const utils::VehicleState& vehicle_state, TrajectoryAnalyzer& ref_trajectory);
    MatrixXd solve_LQR(double tolerance = 0.01, size_t max_iter = 150);
};

double LatController::compute_input(
        const utils::VehicleState& vehicle_state, TrajectoryAnalyzer& ref_trajectory)
{
    Vector4d traj_vec = ref_trajectory.to_trajectory_frame(vehicle_state);
    double theta_e = traj_vec[0];
    double e_cg = traj_vec[1];
    double yaw_ref = traj_vec[2];
    double k_ref = traj_vec[3];

    A(1, 2) = vehicle_state.v;
    B(3, 0) = vehicle_state.v / vehicle_state.vc.WB;
    MatrixXd K = solve_LQR();
    Matrix<double, 4, 1> x = Matrix<double, 4, 1>::Zero();
    x << e_cg, (e_cg - e_cg_old) / dt, theta_e, (theta_e - theta_e_old) / dt;

    MatrixXd ustar = -K * x;
    double steer_angle_feedback = utils::pi_2_pi(ustar(0, 0));
    double steer_angle_feedforward = atan2(vehicle_state.vc.WB * k_ref, 1);
    double steer_angle = steer_angle_feedback + steer_angle_feedforward;
    
    e_cg_old = e_cg;
    theta_e_old = theta_e;

    return steer_angle;
}

MatrixXd LatController::solve_LQR(double tolerance, size_t max_iter)
{
    MatrixXd P = Q;
    MatrixXd P_next = Q;

    // solve a Discrete-time Algebraic Riccati Equation
    for (size_t i = 0; i < max_iter; ++i) {
        P_next = A.transpose() * P * A - A.transpose() * P * B * 
            (R + B.transpose() * P * B).inverse() * B.transpose() * P * A + Q;

        if ((P_next - P).array().abs().maxCoeff() < tolerance) {
            P = P_next;
            break;
        }
        P = P_next;
    }

    MatrixXd K = ((B.transpose() * P * B + R)).inverse() * (B.transpose() * P * A);

    return K;
}

class LonController
{
private:
    double kp;
public:
    explicit LonController(double _kp = 0.3) : kp(_kp) {}
    LonController() = delete;
    ~LonController() {}

    double compute_input(double target_speed, const utils::VehicleState& vehicle_state, double dist);
};

double LonController::compute_input(double target_speed, const utils::VehicleState& vehicle_state, double dist)
{
    // Longitudinal Controller using PID
    // Currently, the planned path is not converted to the frenet coordinate,
    // and this is not a longitudinal control in the frenet coordinate.
    // Fortunately, if there is a trajectory in the frenet coordinate, it can be easily converted.
    double accel = kp * (target_speed - vehicle_state.v);

    if (dist < 10.0) {
        if (vehicle_state.v > 2.0) {
            accel = -3.0;
        } else if (vehicle_state.v < -2) {
            accel = -1.0;
        }
    }

    return accel;
}

int main(int argc, char** argv)
{
    vector<double> ax = {0.0, 10., 16., 20.0, 14., 4, 8};
    vector<double> ay = {0.0, -4., 2., 4.0, 12., 8, 4};
    Vector2d goal(ax.back(), ay.back());
    vector<vector<double>> traj = CubicSpline2D::calc_spline_course(ax, ay, 0.1);
    double time = 0.0;
    double target_speed = 12.0 / 3.6;

    utils::VehicleConfig vc(0.5);
    vc.MAX_STEER = M_PI_4;
    utils::VehicleState state(vc, 0., 0., 0., 0.);

    LatController lat_controller(DT);
    LonController lon_controller(0.3);
    TrajectoryAnalyzer ref_trajectory(traj[0], traj[1], traj[2], traj[3]);
    
    vector<double> x = {state.x};
    vector<double> y = {state.y};
    vector<double> yaw = {state.yaw};
    vector<double> v = {state.v};

    while (time < MAX_SIM_TIME) {
        double dist = hypot(state.x - goal[0], state.y - goal[1]);
        if (dist <= GOAL_DIS) {
            break;
        }

        double steer = lat_controller.compute_input(state, ref_trajectory);
        double acc = lon_controller.compute_input(target_speed, state, dist);
        state.update(acc, steer, DT);

        time = time + DT;
        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);

        if (show_animation) {
            plt::cla();
            plt::named_plot("course", traj[0], traj[1], "-r");
            plt::named_plot("trajectory", x, y, "-b");
            utils::draw_vehicle({state.x, state.y, state.yaw}, steer, vc);
            
            plt::title("LQR with frenet frame Speed [km/h]: " +
                            std::to_string(round(state.v * 3.6)).substr(0, 4));
            plt::axis("equal");
            plt::legend({{"loc", "upper left"}});
            plt::grid(true);
            plt::pause(0.0001);
        }
    }
    plt::legend({{"loc", "upper left"}});
    plt::show();

    return 0;
}
