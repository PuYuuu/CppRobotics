#include <cmath>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <fmt/core.h>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"
#include "PathPlanning/CurvesGenerator/cubic_spline.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr double MAX_SIM_TIME = 500.0;
constexpr double GOAL_DIS = 0.3;
constexpr double STOP_SPEED = 0.05;
constexpr double DT = 0.1;
constexpr bool show_animation = true;

// How to design a universal customizable state vector
// LQR controller is a difficult problem for me.
class LQRController
{
private:
    MatrixXd A;
    MatrixXd B;
    MatrixXd Q;
    MatrixXd R;
    double pe;
    double pth_e;

    MatrixXd solve_LQR(void);
    MatrixXd solve_dare(double tolerance = 0.01, size_t max_iter = 150);
public:
    LQRController(MatrixXd a, MatrixXd b, MatrixXd q, MatrixXd r) :
        A(a), B(b), Q(q), R(r), pe(0.), pth_e(0.) {}
    ~LQRController() {}

    Vector2d compute_input(const utils::VehicleState& state, Vector4d target, double tv);
};

MatrixXd LQRController::solve_LQR(void)
{
    MatrixXd P = solve_dare();
    // compute the LQR gain
    MatrixXd K = ((B.transpose() * P * B + R)).inverse() * (B.transpose() * P * A); 

    return K;
}

// solve a Discrete-time Algebraic Riccati Equation (DARE)
// x_{k+1} = A * x_{k} + B * u_{k}
// J = sum{ x_{k}.T * Q * x_{k} + u_{k}.T * R * u_{k} }
MatrixXd LQRController::solve_dare(double tolerance, size_t max_iter)
{
    MatrixXd p = Q;
    MatrixXd p_next = Q;

    for (size_t i = 0; i < max_iter; ++i) {
        p_next = A.transpose() * p * A - A.transpose() * p * B * 
            (R + B.transpose() * p * B).inverse() * B.transpose() * p * A + Q;

        if ((p_next - p).array().abs().maxCoeff() < tolerance) {
            break;
        }
        p = p_next;
    }

    return p_next;
}

Vector2d LQRController::compute_input(const utils::VehicleState& state, Vector4d target, double tv)
{
    double dxl = target[0] - state.x;
    double dyl = target[1] - state.y;
    double e = hypot(dxl, dyl);
    double angle = utils::pi_2_pi(target[2] - atan2(dyl, dxl));
    if (angle < 0) {
        e *= -1;
    }

    double v = state.v;
    double th_e = utils::pi_2_pi(state.yaw - target[2]);
    A(1, 2) = v;
    B(3, 0) = v / state.vc.WB;
    MatrixXd K = solve_LQR();
    // state vector x = [e, dot_e, th_e, dot_th_e, delta_v]
    Matrix<double, 5, 1> x = Matrix<double, 5, 1>::Zero();
    x << e, (e - pe) / DT, th_e, (th_e - pth_e) / DT, v - tv;
    
    MatrixXd ustar = -K * x;
    double steer_angle_feedforward = atan2(state.vc.WB * target[3], 1);
    double steer_angle_feedback = utils::pi_2_pi(ustar(0, 0));
    double delta = steer_angle_feedforward + steer_angle_feedback;
    double accel = ustar(1, 0);

    pe = e;
    pth_e = th_e;

    return {accel, delta};
}

vector<double> calc_speed_profile(const vector<double>& cyaw, double target_speed)
{
    int len = cyaw.size();
    vector<double> speed_profile(len, target_speed);
    int direction = 1;

    for (size_t idx = 0; idx < len - 1; ++idx) {
        double dyaw = abs(cyaw[idx + 1] - cyaw[idx]);
        bool sw = (M_PI_4 <= dyaw && dyaw < M_PI_2);

        if (sw) {
            direction *= -1;
        }
        if (direction < 0) {
            speed_profile[idx] = -target_speed;
        } else {
            speed_profile[idx] = target_speed;
        }
        if (sw) {
            speed_profile[idx] = 0;
        }
    }

    for (size_t idx = 0; idx < 30; ++idx) {
        speed_profile[len - 1 - idx] = target_speed / (30 - idx);
        if (speed_profile[len - 1 - idx] <= 1.0 / 3.6) {
            speed_profile[len - 1 - idx] = 1.0 / 3.6;
        }
    }

    return speed_profile;
}

size_t calc_nearest_index(const utils::VehicleState& state,
    const vector<double>& cx, const vector<double>& cy)
{
    size_t nearest_index = -1;
    double min_d = 1e6;
    for (size_t idx = 0; idx < cx.size(); ++idx) {
        double dx = state.x - cx[idx];
        double dy = state.y - cy[idx];
        double d = hypot(dx, dy);
        if (nearest_index == -1 || d <= min_d) {
            nearest_index = idx;
            min_d = d;
        }
    }
    
    return nearest_index;
}

int main(int argc, char** argv)
{
    vector<double> ax = {0.0, 10., 16., 20.0, 14., 4, 8};
    vector<double> ay = {0.0, -4., 2., 4.0, 12., 8, 4};
    Vector2d goal(ax.back(), ay.back());
    vector<vector<double>> traj = CubicSpline2D::calc_spline_course(ax, ay, 0.1);

    double target_speed = 10.0 / 3.6;
    double time = 0.0;
    utils::VehicleConfig vc(0.5);
    vc.MAX_STEER = M_PI_4;
    utils::VehicleState state(vc, 0., 0., 0., 0.);
    vector<double> sp = calc_speed_profile(traj[2], target_speed);

    Matrix<double, 5, 5> A = Matrix<double, 5, 5>::Zero();
    A(0, 0) = 1.0;
    A(0, 1) = DT;
    A(1, 2) = state.v;
    A(2, 2) = 1.0;
    A(2, 3) = DT;
    A(4, 4) = 1.0;
    Matrix<double, 5, 2> B = Matrix<double, 5, 2>::Zero();
    B(3, 0) = state.v / state.vc.WB;
    B(4, 1) = DT;
    Matrix<double, 5, 5> Q = Matrix<double, 5, 5>::Identity();
    Matrix<double, 2, 2> R = Matrix<double, 2, 2>::Identity();
    LQRController lqr(A, B, Q, R);

    vector<double> x = {state.x};
    vector<double> y = {state.y};
    vector<double> yaw = {state.yaw};
    vector<double> v = {state.v};

    while (time < MAX_SIM_TIME) {
        size_t ind = calc_nearest_index(state, traj[0], traj[1]);
        Vector2d control = lqr.compute_input(
            state, {traj[0][ind], traj[1][ind], traj[2][ind], traj[3][ind]}, sp[ind]);
        state.update(control[0], control[1], DT);

        time = time + DT;
        if (hypot(state.x - goal[0], state.y - goal[1]) <= GOAL_DIS) {
            break;
        }
        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);

        if (show_animation) {
            plt::cla();
            plt::named_plot("course", traj[0], traj[1], "-k");
            plt::named_plot("trajectory", x, y, "-r");
            plt::plot({traj[0][ind]}, {traj[1][ind]}, "xg");
            utils::draw_vehicle({state.x, state.y, state.yaw}, control[1], vc);
            
            plt::title("LQR with cartesian frame Speed [km/h]: " +
                            std::to_string(round(state.v * 3.6)).substr(0, 4));
            plt::axis("equal");
            plt::legend({{"loc", "upper left"}});
            plt::grid(true);
            plt::pause(0.0001);
        }
    }
    plt::named_plot("waypoint", ax, ay, "xg");
    plt::legend({{"loc", "upper left"}});
    plt::show();

    return 0;
}
