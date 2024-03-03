#include <cmath>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <fmt/core.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

#include "utils.hpp"
#include "matplotlibcpp.h"
#include "PathPlanning/include/cubic_spline.hpp"

using std::vector;
using CppAD::AD;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr size_t NX = 4;  // x = [x, y, v, yaw]
constexpr size_t NU = 2;  // a = [accel, steer]
constexpr size_t TT = 20; // horizon length

constexpr double GOAL_DIS = 1.0;            // goal distance
constexpr double MAX_SIM_TIME = 500.0;      // max simulation time
constexpr double MAX_STEER = M_PI_4;        // maximum steering angle [rad]
constexpr double MAX_DSTEER = M_PI_2 / 3;   // maximum steering speed [rad/s]
constexpr double MAX_SPEED = 55.0 / 3.6;    // maximum speed [m/s]
constexpr double MIN_SPEED = -20.0 / 3.6;   // minimum speed [m/s]
constexpr double MAX_ACCEL = 1.0;           // maximum accel [m/ss]

constexpr double TARGET_SPEED = 10.0 / 3.6; // [m/s] target speed

constexpr double DT = 0.2;      // [s] time tick

constexpr size_t x_start = 0;
constexpr size_t y_start = x_start + TT;
constexpr size_t yaw_start = y_start + TT;
constexpr size_t v_start = yaw_start + TT;

constexpr size_t delta_start = v_start + TT;
constexpr size_t a_start = delta_start + TT - 1;

constexpr double WB = 2.5;
constexpr double show_animation = true;

double find_nearest_s(Vector2d p, const CubicSpline2D& sp) 
{
    static double last_min_s = sp.s.front();
    double min_dist = (sp.calc_position(sp.s.front()) - p).norm();
    double min_s = last_min_s;

    for (double s = last_min_s; s < sp.s.back() && s < last_min_s + 5; s += 0.01) {
        double dist = (sp.calc_position(s) - p).norm();
        if (dist < min_dist) {
            min_s = s;
            min_dist = dist;
        }
    }
    last_min_s = min_s;

    return min_s;
}

void cal_ref_point(double s0, Vector4d& state, const CubicSpline2D& sp) 
{
    Vector2d xy = sp(s0, 0);
    Vector2d dxy = sp(s0, 1);
    Vector2d ddxy = sp(s0, 2);
    double dx = dxy.x();
    double dy = dxy.y();
    double ddx = ddxy.x();
    double ddy = ddxy.y();
    double dphi = (ddy * dx - dy * ddx) / (dx * dx + dy * dy);
    state[0] = xy.x();
    state[1] = xy.y();
    state[2] = atan2(dy, dx);
    state[3] = atan2(WB * dphi, 1.0);
}

MatrixXd calc_ref_trajectory(const utils::VehicleState& state, const CubicSpline2D& sp)
{
    MatrixXd ref_traj = Matrix<double, 5, TT>::Zero();
    double s0 = find_nearest_s({state.x, state.y}, sp);

    for (int i = 0; i < TT; ++i) {
        Vector4d ref_state;
        cal_ref_point(s0, ref_state, sp);
        ref_traj(0, i) = ref_state[0];
        ref_traj(1, i) = ref_state[1];
        ref_traj(2, i) = ref_state[2];
        ref_traj(3, i) = ref_state[3];
        if (sp.s.back() - s0 < TT * TARGET_SPEED * DT) {
            ref_traj(4, i) = 0;
        } else {
            ref_traj(4, i) = TARGET_SPEED;
        }

        s0 += TARGET_SPEED * DT;
        s0 = s0 < sp.s.back() ? s0 : sp.s.back();
    }

    return ref_traj;
}

class FG_EVAL{
public:
    Matrix<double, 5, TT> traj_ref;

    FG_EVAL(Matrix<double, 5, TT> traj_ref){
        this->traj_ref = traj_ref;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector& fg, const ADvector& vars) {
        fg[0] = 0;

        for(size_t idx = 0; idx < TT - 1; ++idx) {
            fg[0] += 0.01 * CppAD::pow(vars[a_start + idx], 2);
            fg[0] += 0.01 * CppAD::pow(vars[delta_start + idx], 2);
        }

        for(size_t idx = 0; idx < TT - 2; ++idx){
            fg[0] += 0.01 * CppAD::pow(vars[a_start + idx + 1] - vars[a_start + idx], 2);
            fg[0] += 1 * CppAD::pow(vars[delta_start + idx + 1] - vars[delta_start + idx], 2);
        }

        // fix the initial state as a constraint
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + yaw_start] = vars[yaw_start];
        fg[1 + v_start] = vars[v_start];

        fg[0] += CppAD::pow(traj_ref(0, 0) - vars[x_start], 2);
        fg[0] += CppAD::pow(traj_ref(1, 0) - vars[y_start], 2);
        // fg[0] += 0.5 * CppAD::pow(traj_ref(2, 0) - vars[yaw_start], 2);
        fg[0] += CppAD::pow(traj_ref(4, 0) - vars[v_start], 2);

        // The rest of the constraints
        for (size_t idx = 0; idx < TT - 1; ++idx) {
            // The state at time t+1 .
            AD<double> x1 = vars[x_start + idx + 1];
            AD<double> y1 = vars[y_start + idx + 1];
            AD<double> yaw1 = vars[yaw_start + idx + 1];
            AD<double> v1 = vars[v_start + idx + 1];

            // The state at time t.
            AD<double> x0 = vars[x_start + idx];
            AD<double> y0 = vars[y_start + idx];
            AD<double> yaw0 = vars[yaw_start + idx];
            AD<double> v0 = vars[v_start + idx];

            // Only consider the actuation at time t.
            AD<double> delta0 = vars[delta_start + idx];
            AD<double> a0 = vars[a_start + idx];

            // constraint with the dynamic model
            fg[2 + x_start + idx] = x1 - (x0 + v0 * CppAD::cos(yaw0) * DT);
            fg[2 + y_start + idx] = y1 - (y0 + v0 * CppAD::sin(yaw0) * DT);
            fg[2 + yaw_start + idx] = yaw1 - (yaw0 + v0 * CppAD::tan(delta0) / WB * DT);
            fg[2 + v_start + idx] = v1 - (v0 + a0 * DT);
            // cost with the ref traj
            fg[0] += CppAD::pow(traj_ref(0, idx + 1) - (x0 + v0 * CppAD::cos(yaw0) * DT), 2);
            fg[0] += CppAD::pow(traj_ref(1, idx + 1) - (y0 + v0 * CppAD::sin(yaw0) * DT), 2);
            // fg[0] += 0.5 * CppAD::pow(traj_ref(2, idx + 1) - (yaw0 + v0 * CppAD::tan(delta0) / WB * DT), 2);
            fg[0] += CppAD::pow(traj_ref(4, 0) - (v0 + a0 * DT), 2);
        }
    }
};

class MPCController
{
private:
    size_t horizon_length;
    size_t n_vars;
    size_t n_constraints;
    std::string optimize_options;
public:
    using Dvector = CPPAD_TESTVECTOR(double);

    MPCController() = delete;
    MPCController(size_t h_len) : horizon_length(h_len) {
        n_vars = horizon_length * 4 + (horizon_length - 1) * 2;
        n_constraints = horizon_length * 4;

        optimize_options = "";
        optimize_options += "Integer print_level     0\n";
        // optimize_options += "Sparse  true         forward\n";
        optimize_options += "Sparse  true            reverse\n";
        optimize_options += "Integer max_iter        600\n";
        optimize_options += "Numeric tol             1e-8\n";
        optimize_options += "Numeric max_cpu_time    1.0\n";
    }
    ~MPCController() {}

    vector<double> mpc_solve(utils::VehicleState& x0, MatrixXd traj_ref);
    Vector2d compute_input(utils::VehicleState& x0, MatrixXd traj_ref);
};

vector<double> MPCController::mpc_solve(utils::VehicleState& x0, MatrixXd traj_ref)
{
    double x = x0.x;
    double y = x0.y;
    double yaw = x0.yaw;
    double v = x0.v;

    Dvector vars(n_vars);
    for (size_t idx = 0; idx < n_vars; ++idx){
        vars[idx] = 0.0;
    }
    vars[x_start] = x;
    vars[y_start] = y;
    vars[yaw_start] = yaw;
    vars[v_start] = v;

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    for (size_t idx = 0; idx < n_vars; ++idx) {
        vars_lowerbound[idx] = -1e7;
        vars_upperbound[idx] = 1e7;
    }
    for (size_t idx = delta_start; idx < delta_start + TT - 1; ++idx) {
        vars_lowerbound[idx] = -MAX_STEER;
        vars_upperbound[idx] = MAX_STEER;
    }
    for (size_t idx = a_start; idx < a_start + TT - 1; ++idx) {
        vars_lowerbound[idx] = -MAX_ACCEL;
        vars_upperbound[idx] = MAX_ACCEL;
    }
    for (size_t idx = v_start; idx < v_start + TT; ++idx) {
        vars_lowerbound[idx] = MIN_SPEED;
        vars_upperbound[idx] = MAX_SPEED;
    }

    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t idx = 0; idx < n_constraints; ++idx) {
        constraints_lowerbound[idx] = 0;
        constraints_upperbound[idx] = 0;
    }
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[yaw_start] = yaw;
    constraints_lowerbound[v_start] = v;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[yaw_start] = yaw;
    constraints_upperbound[v_start] = v;

    FG_EVAL fg_eval(traj_ref);
    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_EVAL>(
      optimize_options, vars, vars_lowerbound, vars_upperbound,
      constraints_lowerbound, constraints_upperbound, fg_eval, solution);

    bool ok = true;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    vector<double> result;
    for (size_t idx = 0 ; idx < n_vars; ++idx) {
        result.push_back(static_cast<double>(solution.x[idx]));
    }

    return result;
}

Vector2d MPCController::compute_input(utils::VehicleState& x0, MatrixXd traj_ref)
{
    vector<double> output = mpc_solve(x0, traj_ref);
    Vector2d input(output[a_start], output[delta_start]);

    return input;
}

int main(int argc, char** argv)
{
    // vector<double> ax = {0.0, 60.0, 125.0, 50.0, 75.0, 35.0, -10.0};
    // vector<double> ay = {0.0, 0.0, 50.0, 65.0, 30.0, 50.0, -20.0};
    vector<double> ax = {0.0, 30.0, 6.0, 20.0, 35.0, 10.0, -1.0};
    vector<double> ay = {0.0, 0.0, 20.0, 35.0, 20.0, 30.0, -2.0};
    Vector2d goal(ax.back(), ay.back());
    CubicSpline2D sp(ax, ay);

    std::vector<std::vector<double>> course(2);
    for (double i = sp.s.front(); i < sp.s.back(); i += 0.05) {
        Eigen::Vector2d tmp = sp(i, 0);
        course[0].emplace_back(tmp.x());
        course[1].emplace_back(tmp.y());
    }

    utils::VehicleConfig vc;
    vc.MAX_STEER = MAX_STEER;
    vc.WB = WB;
    utils::VehicleState state(vc, 0, 0, 0, 0);
    MPCController mpc(TT);

    double time = 0.0;

    vector<double> x_h;
    vector<double> y_h;
    vector<double> v_h;
    vector<double> t_h;

    while (MAX_SIM_TIME >= time) {
        MatrixXd ref_traj = calc_ref_trajectory(state, sp);

        vector<double> output = mpc.mpc_solve(state, ref_traj);
        vector<vector<double>> ooxy(2);
        for (size_t idx = 1; idx < TT; ++idx) {
            ooxy[0].push_back(output[x_start + idx]);
            ooxy[1].push_back(output[y_start + idx]);
        }

        state.update(output[a_start], output[delta_start], DT);
        double steer = output[delta_start];

        float dx = state.x - goal[0];
        float dy = state.y - goal[1];
        if (hypot(dx, dy) <= GOAL_DIS) {
            fmt::print("Goal!\n");
            break;
        }

        x_h.push_back(state.x);
        y_h.push_back(state.y);
        v_h.push_back(state.v);
        t_h.push_back(time);
        time = time + DT;

        if (show_animation) {
            plt::cla();

            plt::named_plot("course", course[0], course[1], "-r");
            plt::named_plot("trajectory", x_h, y_h, "-b");
            plt::named_plot("prediction", ooxy[0], ooxy[1], "xg");

            utils::draw_vehicle({state.x, state.y, state.yaw}, steer, state.vc);
            plt::axis("equal");
            plt::grid(true);
            plt::title("MPC Tracking Speed[km/h]:" + std::to_string(state.v * 3.6).substr(0,4));
            plt::legend({{"loc", "upper left"}});
            plt::pause(0.01);
        }
    }
    plt::show();

    return 0;
}
