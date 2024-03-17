#pragma once
#ifndef __MODEL_PREDICTIVE_CONTROL_HPP
#define __MODEL_PREDICTIVE_CONTROL_HPP

#include <OsqpEigen/OsqpEigen.h>

#include <Eigen/Core>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <string>
#include <vector>

#include "utils.hpp"

static constexpr int NX = 4;  // state x y phi v
static constexpr int NU = 2;  // input a delta

typedef Eigen::Matrix<double, NX, NX> MatrixA;
typedef Eigen::Matrix<double, NX, NU> MatrixB;
typedef Eigen::Vector4d VectorG;
typedef Eigen::Vector4d VectorX;
typedef Eigen::Vector2d VectorU;

class MPCController {
private:
    size_t n_vars;
    size_t n_constraints;
    std::string optimize_options;

    bool init_ = false;

    double ll_;
    double dt_;
    double rho_;
    int N_;
    double rhoN_;

    double v_max_, a_max_, delta_max_, ddelta_max_;

    std::vector<VectorX> predictState_;
    std::vector<VectorU> predictInput_;

    MatrixA Ad_;
    MatrixB Bd_;
    VectorG gd_;
    // x_{k+1} = Ad * x_{k} + Bd * u_k + gd

    Eigen::SparseMatrix<double> P_, q_, A_, l_, u_;
    Eigen::SparseMatrix<double> Cx_, lx_, ux_;  // p, v constrains
    Eigen::SparseMatrix<double> Cu_, lu_, uu_;  // a delta vs constrains
    Eigen::SparseMatrix<double> Qx_;

public:
    using Dvector = CPPAD_TESTVECTOR(double);

    MPCController();
    ~MPCController() {}

    int solve_with_ipopt(utils::VehicleState& x0, Eigen::MatrixXd traj_ref);
    int solve_with_osqp(utils::VehicleState& x0, Eigen::MatrixXd traj_ref);
    void linearization(const double& phi, const double& v, const double& delta);
    void getPredictXU(std::vector<VectorX>& state, std::vector<VectorU>& input) {
        state = predictState_;
        input = predictInput_;
    }
};

#endif
