#include <cmath>
#include <queue>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <OsqpEigen/OsqpEigen.h>

#include "utils.hpp"
#include "matplotlibcpp.h"
#include "PathPlanning/include/cubic_spline.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr size_t TT = 20;                   // horizon length
constexpr double GOAL_DIS = 1.0;            // goal distance
constexpr double MAX_SIM_TIME = 500.0;      // max simulation time
constexpr double MAX_STEER = M_PI_4;        // maximum steering angle [rad]
constexpr double MAX_DSTEER = M_PI_2 / 3;   // maximum steering speed [rad/s]
constexpr double MAX_SPEED = 55.0 / 3.6;    // maximum speed [m/s]
constexpr double MIN_SPEED = -20.0 / 3.6;   // minimum speed [m/s]
constexpr double MAX_ACCEL = 2.0;           // maximum accel [m/ss]

constexpr double TARGET_SPEED = 10.0 / 3.6; // [m/s] target speed

constexpr double DT = 0.2;      // [s] time tick

constexpr double WB = 2.5;
constexpr double show_animation = true;

static constexpr int NX = 4;  // state x y phi v
static constexpr int NU = 2;  // input a delta
typedef Eigen::Matrix<double, NX, NX> MatrixA;
typedef Eigen::Matrix<double, NX, NU> MatrixB;
typedef Eigen::Vector4d VectorG;
typedef Eigen::Vector4d VectorX;
typedef Eigen::Vector2d VectorU;

class MPCController
{
private:
    bool init_ = false;

    double ll_;
    double dt_;
    double rho_;
    int N_;
    double rhoN_;

    double v_max_, a_max_, delta_max_, ddelta_max_;

    vector<VectorX> predictState_;
    vector<VectorU> predictInput_;

    MatrixA Ad_;
    MatrixB Bd_;
    VectorG gd_;
    // x_{k+1} = Ad * x_{k} + Bd * u_k + gd

    Eigen::SparseMatrix<double> P_, q_, A_, l_, u_;
    Eigen::SparseMatrix<double> Cx_, lx_, ux_;  // p, v constrains
    Eigen::SparseMatrix<double> Cu_, lu_, uu_;  // a delta vs constrains
    Eigen::SparseMatrix<double> Qx_;
public:

    MPCController(void);
    ~MPCController(void) {}

    int mpc_solve(utils::VehicleState& x0, MatrixXd traj_ref);
    void linearization(const double& phi, const double& v, const double& delta);
    void getPredictXU(vector<VectorX>& state, vector<VectorU>& input) {
        state = predictState_;
        input = predictInput_;
    }
};

MPCController::MPCController(void)
{
    ll_ = WB;
    dt_ = DT;
    rho_ = 1.0;
    N_ = TT;
    rhoN_ = 1.0;
    v_max_ = MAX_SPEED;
    a_max_ = MAX_ACCEL;
    delta_max_ = MAX_STEER;
    ddelta_max_ = MAX_DSTEER;

    Ad_.setIdentity();  // Ad for instance
    Bd_.setZero();
    Bd_(3, 0) = dt_;
    gd_.setZero();
    // set size of sparse matrices
    P_.resize(NU * N_, NU * N_);
    q_.resize(NU * N_, 1);
    Qx_.resize(NX * N_, NX * N_);
    // stage cost
    Qx_.setIdentity();
    for (int i = 1; i < N_; ++i) {
        Qx_.coeffRef(i * NX - 2, i * NX - 2) = rho_;
        Qx_.coeffRef(i * NX - 1, i * NX - 1) = 0;
    }
    Qx_.coeffRef(N_ * NX - 4, N_ * NX - 4) = rhoN_;
    Qx_.coeffRef(N_ * NX - 3, N_ * NX - 3) = rhoN_;
    Qx_.coeffRef(N_ * NX - 2, N_ * NX - 2) = rhoN_ * rho_;
    int n_cons = 4;  // v a delta ddelta
    A_.resize(n_cons * N_, NU * N_);
    l_.resize(n_cons * N_, 1);
    u_.resize(n_cons * N_, 1);
    // v constrains
    Cx_.resize(1 * N_, NX * N_);
    lx_.resize(1 * N_, 1);
    ux_.resize(1 * N_, 1);
    // a delta constrains
    Cu_.resize(3 * N_, NU * N_);
    lu_.resize(3 * N_, 1);
    uu_.resize(3 * N_, 1);
    // set lower and upper boundaries
    for (int i = 0; i < N_; ++i) {
        // set stage constraints of inputs (a, delta, ddelta)
        // -a_max <= a <= a_max for instance:
        Cu_.coeffRef(i * 3 + 0, i * NU + 0) = 1;
        Cu_.coeffRef(i * 3 + 1, i * NU + 1) = 1;
        Cu_.coeffRef(i * 3 + 2, i * NU + 1) = 1;
        lu_.coeffRef(i * 3 + 0, 0) = -a_max_;
        uu_.coeffRef(i * 3 + 0, 0) = a_max_;
        lu_.coeffRef(i * 3 + 1, 0) = -delta_max_;
        uu_.coeffRef(i * 3 + 1, 0) = delta_max_;
        lu_.coeffRef(i * 3 + 2, 0) = -ddelta_max_ * dt_;
        uu_.coeffRef(i * 3 + 2, 0) = ddelta_max_ * dt_;
        if (i > 0) {
            Cu_.coeffRef(i * 3 + 2, (i - 1) * NU + 1) = -1;
        }

        // set stage constraints of states (v)
        // -v_max <= v <= v_max
        Cx_.coeffRef(i, i * NX + 3) = 1;
        lx_.coeffRef(i, 0) = -0.1;
        ux_.coeffRef(i, 0) = v_max_;
    }
    // set predict mats size
    predictState_.resize(N_);
    predictInput_.resize(N_);
    for (int i = 0; i < N_; ++i) {
        predictInput_[i].setZero();
    }
}

void MPCController::linearization(const double& phi, const double& v, const double& delta)
{
    // set values to Ad_, Bd_, gd_
    Ad_(0, 2) = -v * sin(phi) * dt_;
    Ad_(0, 3) = cos(phi) * dt_;
    Ad_(1, 2) = v * cos(phi) * dt_;
    Ad_(1, 3) = sin(phi) * dt_;
    Ad_(2, 3) = tan(delta) / ll_ * dt_;
    Bd_(2, 1) = v / ll_ / cos(delta) / cos(delta) * dt_;
    gd_(0) = v * sin(phi) * dt_ * phi;
    gd_(1) = -v * cos(phi) * dt_ * phi;
    gd_(2) = -v / ll_ / cos(delta) / cos(delta) * dt_ * delta;
}

int MPCController::mpc_solve(utils::VehicleState& x0_, MatrixXd traj_ref)
{
    lu_.coeffRef(2, 0) = predictInput_.front()(1) - ddelta_max_ * dt_;
    uu_.coeffRef(2, 0) = predictInput_.front()(1) + ddelta_max_ * dt_;
    VectorX x0 = {x0_.x, x0_.y, x0_.yaw, x0_.v};

    MatrixXd BB, AA, gg;
    BB.setZero(NX * N_, NU * N_);
    AA.setZero(NX * N_, NX);
    gg.setZero(NX * N_, 1);
    double phi, v, delta;
    double last_phi = x0(2);
    Eigen::SparseMatrix<double> qx;
    qx.resize(NX * N_, 1);

    for (int i = 0; i < N_; ++i) {
        phi = traj_ref(2, i);
        delta = traj_ref(3, i);
        v = traj_ref(4, i);
        if (phi - last_phi > M_PI) {
            phi -= 2 * M_PI;
        } else if (phi - last_phi < -M_PI) {
            phi += 2 * M_PI;
        }
        last_phi = phi;
        if (init_) {
            double phii = predictState_[i](2);
            v = predictState_[i](3);
            delta = predictInput_[i](1);
            if (phii - last_phi > M_PI) {
                phii -= 2 * M_PI;
            } else if (phii - last_phi < -M_PI) {
                phii += 2 * M_PI;
            }
            last_phi = phii;
            linearization(phii, v, delta);
        } else {
            linearization(phi, v, delta);
        }
        // calculate big state-space matrices
        /* *                BB                AA
        * x1    /       B    0  ... 0 \    /   A \
        * x2    |      AB    B  ... 0 |    |  A2 |
        * x3  = |    A^2B   AB  ... 0 |u + | ... |x0 + gg
        * ...   |     ...  ...  ... 0 |    | ... |
        * xN    \A^(n-1)B  ...  ... B /    \ A^N /
        *
        *     X = BB * U + AA * x0 + gg
        * */
        if (i == 0) {
            BB.block(0, 0, NX, NU) = Bd_;
            AA.block(0, 0, NX, NX) = Ad_;
            gg.block(0, 0, NX, 1) = gd_;
        } else {
            // set BB AA gg
            BB.block(NX * i, 0, NX, NU * N_) = Ad_ * BB.block(NX * (i - 1), 0, NX, NU * N_);
            BB.block(NX * i, NU * i, NX, NU) = Bd_;
            AA.block(NX * i, 0, NX, NX) = Ad_ * AA.block(NX * (i - 1), 0, NX, NX);
            gg.block(NX * i, 0, NX, 1) = Ad_ * gg.block(NX * (i - 1), 0, NX, 1) + gd_;
        }
        // set qx
        qx.coeffRef(i * NX + 0, 0) = -traj_ref(0, i);
        qx.coeffRef(i * NX + 1, 0) = -traj_ref(1, i);
        qx.coeffRef(i * NX + 2, 0) = -rho_ * phi;
        if (i == N_ - 1) {
            qx.coeffRef(i * NX + 0, 0) *= rhoN_;
            qx.coeffRef(i * NX + 1, 0) *= rhoN_;
            qx.coeffRef(i * NX + 2, 0) *= rhoN_;
        }
    }

    Eigen::SparseMatrix<double> BB_sparse = BB.sparseView();
    Eigen::SparseMatrix<double> AA_sparse = AA.sparseView();
    Eigen::SparseMatrix<double> gg_sparse = gg.sparseView();
    Eigen::SparseMatrix<double> x0_sparse = x0.sparseView();

    Eigen::SparseMatrix<double> Cx = Cx_ * BB_sparse;
    Eigen::SparseMatrix<double> lx = lx_ - Cx_ * AA_sparse * x0_sparse - Cx_ * gg_sparse;
    Eigen::SparseMatrix<double> ux = ux_ - Cx_ * AA_sparse * x0_sparse - Cx_ * gg_sparse;
    Eigen::SparseMatrix<double> A_T = A_.transpose();
    A_T.middleCols(0, Cx.rows()) = Cx.transpose();
    A_T.middleCols(Cx.rows(), Cu_.rows()) = Cu_.transpose();
    A_ = A_T.transpose();
    for (int i = 0; i < lx.rows(); ++i) {
        l_.coeffRef(i, 0) = lx.coeff(i, 0);
        u_.coeffRef(i, 0) = ux.coeff(i, 0);
    }
    for (int i = 0; i < lu_.rows(); ++i) {
        l_.coeffRef(i + lx.rows(), 0) = lu_.coeff(i, 0);
        u_.coeffRef(i + lx.rows(), 0) = uu_.coeff(i, 0);
    }
    Eigen::SparseMatrix<double> BBT_sparse = BB_sparse.transpose();
    P_ = BBT_sparse * Qx_ * BB_sparse;
    q_ = BBT_sparse * Qx_.transpose() * (AA_sparse * x0_sparse + gg_sparse) + BBT_sparse * qx;

    // osqp
    OsqpEigen::Solver solver;
    Eigen::VectorXd q_d = q_.toDense();
    Eigen::VectorXd l_d = l_.toDense();
    Eigen::VectorXd u_d = u_.toDense();

    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(P_.rows());
    solver.data()->setNumberOfConstraints(A_.rows());
    solver.data()->setHessianMatrix(P_);
    solver.data()->setGradient(q_d);
    solver.data()->setLinearConstraintsMatrix(A_);
    solver.data()->setLowerBound(l_d);
    solver.data()->setUpperBound(u_d);

    if (!solver.initSolver()) {
        return -1;
    }

    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        return -1;
    }

    VectorXd sol = solver.getSolution();
    MatrixXd solMat = Eigen::Map<const MatrixXd>(sol.data(), NU, N_);

    VectorXd solState = BB * sol + AA * x0 + gg;
    MatrixXd predictMat = Eigen::Map<const MatrixXd>(solState.data(), NX, N_);

    for (int i = 0; i < N_; ++i) {
        predictInput_[i] = solMat.col(i);
        predictState_[i] = predictMat.col(i);
    }
    init_ = true;

    return 0;
}

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

int main(int argc, char** argv)
{
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
    double time = 0.0;

    vector<double> x_h;
    vector<double> y_h;
    vector<double> v_h;
    vector<double> t_h;

    MPCController mpc;
    vector<VectorX> x;
    vector<VectorU> u;
    while (MAX_SIM_TIME >= time) {
        vector<vector<double>> ooxy(2);

        MatrixXd ref_traj = calc_ref_trajectory(state, sp);
        int ret = mpc.mpc_solve(state, ref_traj);
        if (ret != 0) {
            fmt::print("mpc solve error !\n");
        }

        mpc.getPredictXU(x, u);
        for (size_t idx = 1; idx < TT; ++idx) {
            ooxy[0].push_back(x[idx][0]);
            ooxy[1].push_back(x[idx][1]);
        }

        state.update(u.front()[0], u.front()[1], DT);

        double dx = state.x - goal.x();
        double dy = state.y - goal.y();
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

            utils::draw_vehicle({state.x, state.y, state.yaw}, u.front()[1], state.vc);
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
