#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <fmt/core.h>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;
constexpr double DT = 0.1;
constexpr double SIM_TIME = 50.0;
constexpr bool show_animation = true;

Vector2d calc_input(void)
{
    double v = 1.0;         // [m/s]
    double yawrate = 0.1;   // [rad/s]
    
    return {v, yawrate};
}

Vector4d motion_model(Vector4d x, Vector2d u)
{
    Matrix4d F;
    F << 1., 0., 0., 0.,
         0., 1., 0., 0.,
         0., 0., 1., 0.,
         0., 0., 0., 1.;
    
    Matrix<double, 4, 2> B;
    B << DT * cos(x(2,0)),  0,
         DT * sin(x(2,0)),  0,
         0.0,  DT,
         1.0,  0.0;

    Vector4d x_next = F * x + B * u;

    return x_next;
}

Vector2d observation_model(Vector4d x)
{
    Matrix<double, 2, 4> H;
    H << 1, 0, 0, 0,
         0, 1, 0, 0;

    Vector2d z = H * x;

    return z;
}

Matrix4d jacob_f(Vector4d x, Vector2d u)
{
    Matrix4d jF = Matrix4d::Identity();
    double yaw = x(2);
    double v = u(0);
    jF(0, 2) = -DT * v * sin(yaw);
    jF(0, 3) = DT * cos(yaw);
    jF(1, 2) = DT * v * cos(yaw);
    jF(1, 3) = DT * sin(yaw);
    
    return jF;
}

Matrix<double, 2, 4> jacob_h()
{
    Matrix<double, 2, 4> jH;
    jH << 1, 0, 0, 0,
          0, 1, 0, 0;

    return jH;
}

void ekf_estimation(Vector4d& xEst, Matrix4d& PEst,
    Vector2d z, Vector2d u, Matrix4d Q, Matrix2d R)
{
    Vector4d xPred = motion_model(xEst, u);
    Matrix4d jF = jacob_f(xPred, u);
    Matrix4d PPred = jF * PEst * jF.transpose() + Q;

    Matrix<double, 2, 4> jH = jacob_h();
    Vector2d zPred = observation_model(xPred);
    Vector2d y = z - zPred;
    Matrix2d S = jH * PPred * jH.transpose() + R;
    Matrix<double, 4, 2> K = PPred * jH.transpose() * S.inverse();
    
    xEst = xPred + K * y;
    PEst = (Matrix4d::Identity() - K * jH) * PPred;
}

void plot_covariance_ellipse(double x, double y, Matrix2d cov,
                    double chi2 = 3.0, std::string style = "-r")
{
    EigenSolver<Matrix2d> ces(cov);
    Matrix2d e_value = ces.pseudoEigenvalueMatrix();
    Matrix2d e_vector = ces.pseudoEigenvectors();

    double angle = atan2(e_vector(0, 1), e_vector(0, 0));
    Matrix2d l_rot = Utils::rotation_matrix2d(angle);
    vector<double> px;
    vector<double> py;
    double a = e_vector(0, 1);
    double b = e_vector(0, 0);
    if (a < b) {
        std::swap(a, b);
    }

    for (double t = 0; t < 2 * M_PI + 0.2; t += 0.2) {
        Vector2d pxy;
        pxy << 2.0 * a * cos(t), 2.0 * b * sin(t);
        Matrix<double, 1, 2> l_pxy = pxy.transpose() * l_rot;
        px.push_back(l_pxy(0, 0) + x);
        py.push_back(l_pxy(0, 1) + y);
    }
    plt::plot(px, py, style);
}

int main(int argc, char** argv)
{
    double time = 0.0;
    Vector2d u(1., 0.1); // control input
    Vector2d ud;         // nosie control input
    Vector2d z;          // observation z

    // State vector [x y yaw v]
    Vector4d xEst(0., 0., 0., 0.);
    Vector4d xTrue(0., 0., 0., 0.);
    Vector4d xDR(0., 0., 0., 0.);
    Matrix4d PEst = Matrix4d::Identity();
    
    vector<Vector4d> hxEst = {xEst};
    vector<Vector4d> hxTrue = {xTrue};
    vector<Vector4d> hxDR = {xTrue};
    vector<Vector2d> hz = {{0., 0.}};

    // Motional model covariance
    Matrix4d Q = Matrix4d::Identity();
    Q(0, 0) = 0.1 * 0.1;
    Q(1, 1) = 0.1 * 0.1;
    Q(2, 2) = (1.0 / 180 * M_PI) * (1.0 / 180 * M_PI);
    Q(3, 3) = 0.1 * 0.1;

    // Observation model covariance
    Matrix2d R = Matrix2d::Identity();

    // Motion model simulation error
    Matrix2d input_noise = Matrix2d::Identity();
    input_noise(1, 1) = (30.0 / 180 * M_PI) * (30.0 / 180 * M_PI);

    // Observation model simulation error
    Matrix2d gps_noise = Matrix2d::Identity();
    gps_noise(0, 0) = 0.5 * 0.5;
    gps_noise(1, 1) = 0.5 * 0.5;

    std::random_device seed;
    std::mt19937 gen(seed());
    std::normal_distribution<> gaussian_d(0, 1);

    while (SIM_TIME >= time) {
        time += DT;
        ud(0) = u(0) + gaussian_d(gen) * input_noise(0, 0);
        ud(1) = u(1) + gaussian_d(gen) * input_noise(1, 1);

        xTrue = motion_model(xTrue, u);
        xDR = motion_model(xDR, ud);

        z(0) = xTrue(0) + gaussian_d(gen) * gps_noise(0,0);
        z(1) = xTrue(1) + gaussian_d(gen) * gps_noise(1,1);

        ekf_estimation(xEst, PEst, z, ud, Q, R);

        hxDR.push_back(xDR);
        hxTrue.push_back(xTrue);
        hxEst.push_back(xEst);
        hz.push_back(z);

        if (show_animation) {
            plt::cla();
            vector<vector<double>> hxDR_vec(2);
            vector<vector<double>> hxTrue_vec(2);
            vector<vector<double>> hxEst_vec(2);
            vector<vector<double>> hz_vec(2);
            for (size_t idx = 0; idx < hxDR.size(); ++idx) {
                hxDR_vec[0].push_back(hxDR[idx](0));
                hxDR_vec[1].push_back(hxDR[idx](1));
                hxTrue_vec[0].push_back(hxTrue[idx](0));
                hxTrue_vec[1].push_back(hxTrue[idx](1));
                hxEst_vec[0].push_back(hxEst[idx](0));
                hxEst_vec[1].push_back(hxEst[idx](1));
                hz_vec[0].push_back(hz[idx](0));
                hz_vec[1].push_back(hz[idx](1));
            }

            plt::named_plot("GPS", hz_vec[0], hz_vec[1], ".g");
            plt::named_plot("Ground-truth", hxTrue_vec[0], hxTrue_vec[1], "-b");
            plt::named_plot("Dead-reckoning", hxDR_vec[0], hxDR_vec[1], "-k");
            plt::named_plot("Estimation", hxEst_vec[0], hxEst_vec[1], "-r");
            plot_covariance_ellipse(xEst(0), xEst(1), PEst.block(0,0,2,2));
            plt::title("EKF Localization");
            plt::axis("equal");
            plt::grid(true);
            plt::legend({{"loc", "upper right"}});
            plt::pause(0.005);
        }
    }
    plt::show();

    return 0;
}
