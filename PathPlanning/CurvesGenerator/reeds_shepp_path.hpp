#pragma once
#ifndef __REEDS_SHEPP_PATH_HPP
#define __REEDS_SHEPP_PATH_HPP

#include <vector>
#include <Eigen/Core>

class Path
{
public:
    std::vector<double> lengths;
    std::vector<char> ctypes;
    double L;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> yaw;
    std::vector<double> yawt;
    std::vector<int> directions;

    Path(std::vector<double> _x, std::vector<double> _y, std::vector<double> _yaw,
        std::vector<int> _dir) : x(_x), y(_y), yaw(_yaw), directions(_dir) {}
    Path(std::vector<double> _x, std::vector<double> _y, std::vector<double> _yaw,
        std::vector<double> _yawt, std::vector<int> _dir) : x(_x), y(_y), yaw(_yaw),
        yawt(_yawt), directions(_dir) {}
    Path() {}
    ~Path() {}
};

std::vector<Path> calc_rs_paths(
    Eigen::Vector3d s, Eigen::Vector3d g, double maxc, double step_size);
Path reeds_shepp_path(
    Eigen::Vector3d s, Eigen::Vector3d g, double maxc, double step_size = 0.2);

#endif
