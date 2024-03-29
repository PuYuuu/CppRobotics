#pragma once
#ifndef __SIMULATOR_HPP
#define __SIMULATOR_HPP

#include <Eigen/Core>
#include <cmath>
#include <random>
#include <vector>

#include "matplotlibcpp.h"
#include "utils.hpp"

class VehicleSimulator {
private:
    double x;
    double y;
    double yaw;
    double v;
    double max_v;
    double W;
    double L;
    int id;
    std::vector<Eigen::Vector2d> vc_xy;
    static int global_id;

public:
    VehicleSimulator(double _x, double _y, double _yaw, double _v, double _max_v, double _w,
                     double _l)
        : x(_x), y(_y), yaw(_yaw), v(_v), max_v(_max_v), W(_w), L(_l) {
        vc_xy = {{L / 2.0, W / 2.0},
                 {L / 2.0, -W / 2.0},
                 {-L / 2.0, -W / 2.0},
                 {-L / 2.0, W / 2.0},
                 {L / 2.0, W / 2.0}};

        vc_xy = interpolate(vc_xy);
        id = global_id;
        global_id++;
    }
    ~VehicleSimulator() {}

    std::vector<Eigen::Vector2d> interpolate(std::vector<Eigen::Vector2d>& xy);
    void update(double dt, double a, double omega);
    std::vector<std::vector<double>> calc_global_contour(void);
    void plot(void);
};

int VehicleSimulator::global_id = 0;

std::vector<Eigen::Vector2d> VehicleSimulator::interpolate(std::vector<Eigen::Vector2d>& xy_vec) {
    std::vector<Eigen::Vector2d> rxy;
    double d_theta = 0.05;

    for (size_t idx = 0; idx < xy_vec.size() - 1; ++idx) {
        for (double theta = 0; theta <= 1.0; theta += d_theta) {
            Eigen::Vector2d xy;
            xy[0] = (1.0 - theta) * xy_vec[idx][0] + theta * xy_vec[idx + 1][0];
            xy[1] = (1.0 - theta) * xy_vec[idx][1] + theta * xy_vec[idx + 1][1];
            rxy.emplace_back(xy);
        }
    }

    return rxy;
}

void VehicleSimulator::update(double dt, double a, double omega) {
    x += v * cos(yaw) * dt;
    y += v * sin(yaw) * dt;
    yaw += omega * dt;
    v += a * dt;
    if (v >= max_v) {
        v = max_v;
    }
}

std::vector<std::vector<double>> VehicleSimulator::calc_global_contour(void) {
    std::vector<std::vector<double>> gxy(2);

    for (size_t idx = 0; idx < vc_xy.size(); ++idx) {
        Eigen::Vector2d xy = vc_xy[idx];
        Eigen::Matrix2d g_rot = utils::rotation_matrix2d(yaw);
        Eigen::Matrix<double, 1, 2> g_xy = xy.transpose() * g_rot;
        gxy[0].emplace_back(g_xy(0, 0) + x);
        gxy[1].emplace_back(g_xy(0, 1) + y);
    }

    return gxy;
}

void VehicleSimulator::plot(void) {
    matplotlibcpp::plot({x}, {y}, ".b");
    std::vector<std::vector<double>> gxy = calc_global_contour();
    if (id < 1) {
        matplotlibcpp::named_plot("vehicle", gxy[0], gxy[1], "--b");
    } else {
        matplotlibcpp::plot(gxy[0], gxy[1], "--b");
    }
}

class LidarSimulator {
private:
    double range_noise;
    std::mt19937 engine;

public:
    LidarSimulator(void) : range_noise(0.01) {
        std::random_device seed;
        engine.seed(seed());
    }
    ~LidarSimulator() {}

    std::vector<std::vector<double>> get_observation_points(std::vector<VehicleSimulator> v_list,
                                                            double angle_resolution);
    std::vector<std::vector<double>> ray_casting_filter(std::vector<double> theta_l,
                                                        std::vector<double> range_l,
                                                        double angle_resolution);
};

std::vector<std::vector<double>> LidarSimulator::get_observation_points(
    std::vector<VehicleSimulator> v_list, double angle_resolution) {
    std::vector<double> angle;
    std::vector<double> r;
    std::normal_distribution<> gaussian_d(0, range_noise);

    for (VehicleSimulator v : v_list) {
        std::vector<std::vector<double>> gxy = v.calc_global_contour();
        for (size_t idx = 0; idx < gxy[0].size(); ++idx) {
            double v_angle = atan2(gxy[1][idx], gxy[0][idx]);
            double vr = hypot(gxy[1][idx], gxy[0][idx]) * (1 + gaussian_d(engine));
            angle.push_back(v_angle);
            r.push_back(vr);
        }
    }
    std::vector<std::vector<double>> rxy = ray_casting_filter(angle, r, angle_resolution);

    return rxy;
}

std::vector<std::vector<double>> LidarSimulator::ray_casting_filter(std::vector<double> theta_l,
                                                                    std::vector<double> range_l,
                                                                    double angle_resolution) {
    std::vector<std::vector<double>> rxy(2);
    std::vector<double> range_db(static_cast<int>(floor(M_PI * 2 / angle_resolution) + 1),
                                 std::numeric_limits<double>::max());

    for (size_t idx = 0; idx < theta_l.size(); ++idx) {
        int angle_id = static_cast<int>(round(theta_l[idx] / angle_resolution));
        if (range_db[angle_id] > range_l[idx]) {
            range_db[angle_id] = range_l[idx];
        }
    }

    for (size_t idx = 0; idx < range_db.size(); ++idx) {
        double t = idx * angle_resolution;
        if (range_db[idx] < std::numeric_limits<double>::max()) {
            rxy[0].push_back(range_db[idx] * cos(t));
            rxy[1].push_back(range_db[idx] * sin(t));
        }
    }

    return rxy;
}

#endif
