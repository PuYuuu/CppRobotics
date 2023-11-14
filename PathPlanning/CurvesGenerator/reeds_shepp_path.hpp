#pragma once
#ifndef __REEDS_SHEPP_PATH_HPP
#define __REEDS_SHEPP_PATH_HPP

#include <vector>
#include <Eigen/Core>

using std::vector;
using namespace Eigen;

class Path
{
public:
    vector<double> lengths;
    vector<char> ctypes;
    double L;
    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    vector<int> directions;

    Path(vector<double> _x, vector<double> _y, vector<double> _yaw,
        vector<int> _dir) : x(_x), y(_y), yaw(_yaw), directions(_dir) {

    }
    Path() {}
    ~Path() {}
};

vector<Path> calc_rs_paths(Vector3d s, Vector3d g, double maxc, double step_size);
Path reeds_shepp_path(Vector3d s, Vector3d g, double maxc, double step_size = 0.2);

#endif
