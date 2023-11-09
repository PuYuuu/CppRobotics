#include "utils.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

static void plot(VectorXd vec_x, VectorXd vec_y, std::string style)
{
    std::vector<double> x;
    std::vector<double> y;
    double* x_pointer = vec_x.data();
    double* y_pointer = vec_y.data();
    x.assign(x_pointer, x_pointer + vec_x.rows());
    y.assign(y_pointer, y_pointer + vec_y.rows());

    plt::plot(x, y, style);
}

namespace utils {
    
void draw_vehicle(Vector3d state, double steer, VehicleConfig c, bool draw_wheel, std::string color)
{
    Matrix<double, 2, 5> vehicle;
    Matrix<double, 2, 5> wheel;
    vehicle << -c.RB, -c.RB, c.RF, c.RF, -c.RB,
                c.W / 2, -c.W / 2, -c.W / 2, c.W / 2, c.W / 2;
    wheel << -c.TR, -c.TR, c.TR, c.TR, -c.TR,
            c.TW / 4, -c.TW / 4, -c.TW / 4, c.TW / 4, c.TW / 4;

    Matrix<double, 2, 5> rlWheel = wheel;
    Matrix<double, 2, 5> rrWheel = wheel;
    Matrix<double, 2, 5> frWheel = wheel;
    Matrix<double, 2, 5> flWheel = wheel;

    Matrix2d rot1;
    Matrix2d rot2;
    rot1 << cos(state[2]), -sin(state[2]), sin(state[2]), cos(state[2]);
    rot2 << cos(steer), -sin(steer), sin(steer), cos(steer);

    vehicle = rot1 * vehicle;
    vehicle += Vector2d(state[0], state[1]).replicate(1, 5);
    plot(vehicle.row(0), vehicle.row(1), color);

    if (draw_wheel) {
        frWheel = rot2 * frWheel;
        flWheel = rot2 * flWheel;
    
        frWheel += Vector2d(c.WB, -c.WD / 2).replicate(1, 5);
        flWheel += Vector2d(c.WB, c.WD / 2).replicate(1, 5);
    
        rrWheel.row(1) -= VectorXd::Constant(5, c.WD / 2);
        rlWheel.row(1) += VectorXd::Constant(5, c.WD / 2);
    
        frWheel = rot1 * frWheel;
        flWheel = rot1 * flWheel;
        rrWheel = rot1 * rrWheel;
        rlWheel = rot1 * rlWheel;
        
        frWheel += Vector2d(state[0], state[1]).replicate(1, 5);
        flWheel += Vector2d(state[0], state[1]).replicate(1, 5);
        rrWheel += Vector2d(state[0], state[1]).replicate(1, 5);
        rlWheel += Vector2d(state[0], state[1]).replicate(1, 5);
        
        plot(frWheel.row(0), frWheel.row(1), color);
        plot(flWheel.row(0), flWheel.row(1), color);
        plot(rrWheel.row(0), rrWheel.row(1), color);
        plot(rlWheel.row(0), rlWheel.row(1), color);
    }
}

}
