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
constexpr size_t T  = 6;  // horizon length

constexpr double GOAL_DIS = 1.5;            // goal distance
constexpr double MAX_SIM_TIME = 500.0;      // max simulation time
constexpr double MAX_STEER = M_PI_4;        // maximum steering angle [rad]
constexpr double MAX_DSTEER = M_PI_2 / 3;   // maximum steering speed [rad/s]
constexpr double MAX_SPEED = 55.0 / 3.6;    // maximum speed [m/s]
constexpr double MIN_SPEED = -20.0 / 3.6;   // minimum speed [m/s]
constexpr double MAX_ACCEL = 1.0;           // maximum accel [m/ss]

constexpr double TARGET_SPEED = 10.0 / 3.6; // [m/s] target speed
constexpr double N_IND_SEARCH = 10;         // Search index number

constexpr double DT = 0.2;      // [s] time tick

constexpr size_t x_start = 0;
constexpr size_t y_start = x_start + T;
constexpr size_t yaw_start = y_start + T;
constexpr size_t v_start = yaw_start + T;

constexpr size_t delta_start = v_start + T;
constexpr size_t a_start = delta_start + T - 1;

constexpr double WB = 2.5;
constexpr double show_animation = true;

vector<double> calc_speed_profile(
    const vector<double>& cx, const vector<double>& cy,
    const vector<double>& cyaw, double target_speed)
{
    vector<double> speed_profile(cx.size(), target_speed);

    int direction = 1.0;
    for (size_t idx = 0; idx + 1 < cx.size(); ++idx) {
        double dx = cx[idx + 1] - cx[idx];
        double dy = cy[idx + 1] - cy[idx];
        double move_direction = atan2(dy, dx);

        if (dx != 0. && dy != 0.) {
            double dangle = abs(utils::pi_2_pi(move_direction - cyaw[idx]));
            if (dangle > M_PI_4) {
                direction = -1;
            } else {
                direction = 1;
            }
        }
        if (direction != 1) {
            speed_profile[idx] = -target_speed;
        } else {
            speed_profile[idx] = target_speed;
        }
        // speed_profile[idx] = target_speed;
    }
    speed_profile[cx.size() - 1] = 0.0;

    return speed_profile;
}

int calc_nearest_index(const utils::VehicleState& state, const vector<double>& cx,
    const vector<double>& cy, const vector<double>& cyaw, int pind)
{
    double min_d = 0;
    int min_ind = -1;
    for (size_t idx = pind; idx < pind + N_IND_SEARCH; ++idx) {
        double dx = state.x - cx[idx];
        double dy = state.y - cy[idx];
        double d = hypot(dx, dy);

        if (min_ind == -1 || d < min_d) {
            min_ind = idx;
            min_d = d;
        }
    }

    return min_ind;
}

void smooth_yaw(vector<double>& yaw)
{
    for (size_t idx = 0; idx + 1 < yaw.size(); ++idx) {
        double dyaw = yaw[idx + 1] - yaw[idx];

        while (dyaw >= M_PI_2) {
            yaw[idx + 1] -= (2 * M_PI);
            dyaw = yaw[idx + 1] - yaw[idx];
        }
        while (dyaw <= -M_PI_2) {
            yaw[idx + 1] += (2 * M_PI);
            dyaw = yaw[idx + 1] - yaw[idx];
        }
    }
}

MatrixXd calc_ref_trajectory(const utils::VehicleState& state, const vector<double>& cx,
    const vector<double>& cy, const vector<double>& cyaw, const vector<double>& sp, int& pind)
{
    MatrixXd xref = Matrix<double, NX, T>::Zero();
    size_t ncourse = cx.size();

    int ind = calc_nearest_index(state, cx, cy, cyaw, pind);
    if (pind >= ind) {
        ind = pind;
    }

    xref(0, 0) = cx[ind];
    xref(1, 0) = cy[ind];
    xref(2, 0) = sp[ind];
    xref(3, 0) = cyaw[ind];

    double travel = 0.0;

    for (size_t idx = 0; idx < T; ++idx) {
        travel += (abs(state.v) * DT);
        int dind = static_cast<int>(std::round(travel));

        if ((ind + dind) < ncourse) {
            xref(0, idx) = cx[ind + dind];
            xref(1, idx) = cy[ind + dind];
            xref(2, idx) = cyaw[ind + dind];
            xref(3, idx) = sp[ind + dind];
        } else {
            xref(0, idx) = cx[ncourse - 1];
            xref(1, idx) = cy[ncourse - 1];
            xref(2, idx) = cyaw[ncourse - 1];
            xref(3, idx) = sp[ncourse - 1];
        }
    }
    pind = ind;

    return xref;
}

class FG_EVAL{
public:
    Matrix<double, NX, T> traj_ref;

    FG_EVAL(Matrix<double, NX, T> traj_ref){
        this->traj_ref = traj_ref;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector& fg, const ADvector& vars) {
        fg[0] = 0;

        for(size_t idx = 0; idx < T - 1; ++idx) {
            fg[0] +=  0.01 * CppAD::pow(vars[a_start + idx], 2);
            fg[0] += 0.01 * CppAD::pow(vars[delta_start + idx], 2);
        }

        for(size_t idx = 0; idx < T - 2; ++idx){
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
        fg[0] += 0.5 * CppAD::pow(traj_ref(2, 0) - vars[yaw_start], 2);
        fg[0] += 0.5 * CppAD::pow(traj_ref(3, 0) - vars[v_start], 2);

        // The rest of the constraints
        for (size_t idx = 0; idx < T - 1; ++idx) {
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
            fg[0] += 0.5 * CppAD::pow(traj_ref(2, idx + 1) - (yaw0 + v0 * CppAD::tan(delta0) / WB * DT), 2);
            fg[0] += 0.5 * CppAD::pow(traj_ref(3, idx + 1) - (v0 + a0 * DT), 2);
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
        optimize_options += "Integer max_iter        50\n";
        // optimize_options += "Numeric tol          1e-6\n";
        optimize_options += "Numeric max_cpu_time    0.05\n";
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
    for (size_t idx = delta_start; idx < delta_start + T - 1; ++idx) {
        vars_lowerbound[idx] = -MAX_STEER;
        vars_upperbound[idx] = MAX_STEER;
    }
    for (size_t idx = a_start; idx < a_start + T - 1; ++idx) {
        vars_lowerbound[idx] = -MAX_ACCEL;
        vars_upperbound[idx] = MAX_ACCEL;
    }
    for (size_t idx = v_start; idx < v_start + T; ++idx) {
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
    vector<vector<double>> traj = CubicSpline2D::calc_spline_course(ax, ay, 1.0);
    vector<double> sp = calc_speed_profile(traj[0], traj[1], traj[2], TARGET_SPEED);

    utils::VehicleConfig vc;
    vc.MAX_STEER = MAX_STEER;
    vc.WB = WB;
    utils::VehicleState state(vc, traj[0][0], traj[1][0], traj[2][0], sp[0]);
    MPCController mpc(T);

    if ((state.yaw - traj[2][0]) >= M_PI) {
         state.yaw -= M_PI * 2.0; 
    } else if ((state.yaw - traj[2][0]) <= -1.0 * M_PI) {
        state.yaw += M_PI * 2.0;
    }

    double time = 0.0;
    vector<double> x{state.x};
    vector<double> y{state.y};
    vector<double> yaw{state.yaw};
    vector<double> v{state.v};
    vector<double> t{0.0};
    vector<double> d{0.0};
    vector<double> a{0.0};

    int target_ind = calc_nearest_index(state, traj[0], traj[1], traj[2], 0);
    smooth_yaw(traj[2]);

    vector<double> x_h;
    vector<double> y_h;
    vector<double> v_h;
    vector<double> t_h;

    while (MAX_SIM_TIME >= time) {
        MatrixXd xref =
            calc_ref_trajectory(state, traj[0], traj[1], traj[2], sp, target_ind);
        vector<double> output = mpc.mpc_solve(state, xref);
        vector<vector<double>> ooxy(2);
        for (size_t idx = 1; idx < T; ++idx) {
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
            
            plt::named_plot("course", traj[0], traj[1], "-r");
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
