#include <cmath>

#include "utils/utils.hpp"
#include "reeds_shepp_path.hpp"

using std::vector;
using namespace Eigen;

Vector2d polar(double x, double y)
{
    double r = hypot(x, y);
    double theta = atan2(y, x);

    return {r, theta};
}

Vector3d straight_left_straight(double x, double y, double phi, bool& flag)
{
    flag = false;
    phi = utils::pi_2_pi(phi);
    if (M_PI * 0.01 < phi && phi < M_PI * 0.99 && y != 0) {
        double xd = - y / tan(phi) + x;
        double t = xd - tan(phi / 2.0);
        double u = phi;
        double v = utils::sign(y) * hypot(x - xd, y) - tan(phi / 2.0);
        flag = true;
        return {t, u, v};
    }

    return {0.0, 0.0, 0.0};
}

Vector3d left_straight_left(double x, double y, double phi, bool& flag)
{
    flag = false;
    Vector2d ut = polar(x - sin(phi), y - 1.0 + cos(phi));
    if (ut(1) >= 0.0) {
        double v = utils::pi_2_pi(phi - ut(1));
        if (v >= 0.0) {
            flag = true;
            return {ut(1), ut(0), v};
        }
    }

    return {0.0, 0.0, 0.0};
}

Vector3d left_straight_right(double x, double y, double phi, bool& flag)
{
    flag = false;
    Vector2d ut1 = polar(x + sin(phi), y - 1.0 - cos(phi));
    double u1 = ut1(0) * ut1(0);
    if (u1 >= 4.0) {
        double u = sqrt(u1 - 4.0);
        double theta = atan2(2.0, u);
        double t = utils::pi_2_pi(ut1(1) + theta);
        double v = utils::pi_2_pi(t - phi);

        if (t >= 0.0 && v >= 0.0) {
            flag = true;
            return {t, u, v};
        }
    }

    return {0.0, 0.0, 0.0};
}

Vector3d left_right_left(double x, double y, double phi, bool& flag)
{
    flag = false;
    Vector2d ut1 = polar(x - sin(phi), y - 1.0 + cos(phi));

    if (ut1(0) <= 4.0) {
        double u = -2.0 * asin(0.25 * ut1(0));
        double t = utils::pi_2_pi(ut1(1) + 0.5 * u + M_PI);
        double v = utils::pi_2_pi(phi - t + u);

        if (t >= 0.0 && 0.0 >= u) {
            flag = true;
            return {t, u, v};
        }
    }

    return {0.0, 0.0, 0.0};
}

void set_path(vector<Path>& paths, Vector3d lengths, vector<char> ctypes, double step_size)
{
    Path path;
    path.ctypes = ctypes;
    path.lengths = {lengths(0), lengths(1), lengths(2)};
    path.L = abs(lengths(0)) + abs(lengths(1)) + abs(lengths(2));

    for (Path i_path : paths) {
        bool type_is_same = (i_path.ctypes == path.ctypes);
        bool length_is_close = ((i_path.L - path.L) <= step_size);
        if (type_is_same && length_is_close) {
            return ;
        }
    }

    if (path.L <= step_size) {
        return ;
    }

    paths.push_back(path);
}

void straight_curve_straight(double x, double y, double phi, vector<Path>& paths, double step_size)
{
    bool flag = false;
    Vector3d tuv = straight_left_straight(x, y, phi, flag);
    if (flag) {
        set_path(paths, tuv, {'S', 'L', 'S'}, step_size);
    }

    tuv = straight_left_straight(x, -y, -phi, flag);
    if (flag) {
        set_path(paths, tuv, {'S', 'R', 'S'}, step_size);
    }
}

void curve_straight_curve(double x, double y, double phi, vector<Path>& paths, double step_size)
{
    bool flag = false;
    Vector3d tuv = left_straight_left(x, y, phi, flag);
    if (flag) {
        set_path(paths, tuv, {'L', 'S', 'L'}, step_size);
    }

    tuv = left_straight_left(-x, y, -phi, flag);
    if (flag) {
        set_path(paths, -1 * tuv, {'L', 'S', 'L'}, step_size);
    }

    tuv = left_straight_left(x, -y, -phi, flag);
    if (flag) {
        set_path(paths, tuv, {'R', 'S', 'R'}, step_size);
    }
    tuv = left_straight_left(-x, -y, phi, flag);
    if (flag) {
        set_path(paths, -1 * tuv, {'R', 'S', 'R'}, step_size);
    }

    tuv = left_straight_right(x, y, phi, flag);
    if (flag) {
        set_path(paths, tuv, {'L', 'S', 'R'}, step_size);
    }

    tuv = left_straight_right(-x, y, -phi, flag);
    if (flag) {
        set_path(paths, -1 * tuv, {'L', 'S', 'R'}, step_size);
    }

    tuv = left_straight_right(x, -y, -phi, flag);
    if (flag) {
        set_path(paths, tuv, {'R', 'S', 'L'}, step_size);
    }

    tuv = left_straight_right(-x, -y, phi, flag);
    if (flag) {
        set_path(paths, -1 * tuv, {'R', 'S', 'L'}, step_size);
    }
}

void curve_curve_curve(double x, double y, double phi, vector<Path>& paths, double step_size)
{
    bool flag = false;
    Vector3d tuv = left_right_left(x, y, phi, flag);
    if (flag) {
        set_path(paths, tuv, {'L', 'R', 'L'}, step_size);
    }

    tuv = left_right_left(-x, y, -phi, flag);
    if (flag) {
        set_path(paths, -1 * tuv, {'L', 'R', 'L'}, step_size);
    }

    tuv = left_right_left(x, -y, -phi, flag);
    if (flag) {
        set_path(paths, tuv, {'R', 'L', 'R'}, step_size);
    }

    tuv = left_right_left(-x, -y, phi, flag);
    if (flag) {
        set_path(paths, -1 * tuv, {'R', 'L', 'R'}, step_size);
    }

    double xb = x * cos(phi) + y * sin(phi);
    double yb = x * sin(phi) - y * cos(phi);

    tuv = left_right_left(xb, yb, phi, flag);
    if (flag) {
        set_path(paths, {tuv(2), tuv(1), tuv(0)}, {'L', 'R', 'L'}, step_size);
    }

    tuv = left_right_left(-xb, yb, -phi, flag);
    if (flag) {
        set_path(paths, {-tuv(2), -tuv(1), -tuv(0)}, {'L', 'R', 'L'}, step_size);
    }

    tuv = left_right_left(xb, -yb, -phi, flag);
    if (flag) {
        set_path(paths, {tuv(2), tuv(1), tuv(0)}, {'R', 'L', 'R'}, step_size);
    }

    tuv = left_right_left(-xb, -yb, phi, flag);
    if (flag) {
        set_path(paths, {-tuv(2), -tuv(1), -tuv(0)}, {'R', 'L', 'R'}, step_size);
    }
}

vector<Path> generate_path(Vector3d q0, Vector3d q1, double max_curvature, double step_size)
{
    double dx = q1(0) - q0(0);
    double dy = q1(1) - q0(1);
    double dth = q1(2) - q0(2);
    double c = cos(q0(2));
    double s = sin(q0(2));
    double x = (c * dx + s * dy) * max_curvature;
    double y = (-s * dx + c * dy) * max_curvature;

    vector<Path> paths;
    straight_curve_straight(x, y, dth, paths, step_size);
    curve_straight_curve(x, y, dth, paths, step_size);
    curve_curve_curve(x, y, dth, paths, step_size);

    return paths;
}

vector<vector<double>> calc_interpolate_dists_list(vector<double> lengths, double step_size)
{
    vector<vector<double>> interpolate_dists_list;
    int idx = 0;
    for (double length : lengths) {
        vector<double> interp_dists;
        int len_sign = utils::sign(length);
        for (double d = 0; d < abs(length); d += step_size) {
            interp_dists.push_back(len_sign * d);
        }
        interp_dists.push_back(length);
        interpolate_dists_list.emplace_back(interp_dists);
    }

    return interpolate_dists_list;
}

Vector4d interpolate(double dist, double length, char mode, double max_curvature, Vector3d origin)
{
    Vector4d inter(0, 0, 0, 0);
    if (mode == 'S') {
        inter(0) = origin(0) + dist / max_curvature * cos(origin(2));
        inter(1) = origin(1) + dist / max_curvature * sin(origin(2));
        inter(2) = origin(2);
    } else {
        double ldx = sin(dist) / max_curvature;
        double ldy = 0.0;
        if (mode == 'L') {
            ldy = (1.0 - cos(dist)) / max_curvature;
            inter(2) = origin(2) + dist;
        } else if (mode == 'R') {
            ldy = (1.0 - cos(dist)) / -max_curvature;
            inter(2) = origin(2) - dist;
        }
        double gdx = cos(-origin(2)) * ldx + sin(-origin(2)) * ldy;
        double gdy = -sin(-origin(2)) * ldx + cos(-origin(2)) * ldy;
        inter(0) = origin(0) + gdx;
        inter(1) = origin(1) + gdy;
    }
    inter(3) = length > 0.0 ? 1 : -1;
    
    return inter;
}

vector<vector<double>> generate_local_course(
    vector<double> lengths, vector<char> modes, double max_curvature, double step_size)
{
    
    vector<vector<double>> interpolate_dists_list = 
                            calc_interpolate_dists_list(lengths, step_size);

    Vector3d origin(0, 0, 0);
    vector<double> xs;
    vector<double> ys;
    vector<double> yaws;
    vector<double> directions;
    
    for (size_t idx = 0; idx < lengths.size(); ++idx) {
        for (double dist : interpolate_dists_list[idx]) {
            Vector4d state = interpolate(dist, lengths[idx], modes[idx],
                                        max_curvature, origin);
            xs.push_back(state(0));
            ys.push_back(state(1));
            yaws.push_back(state(2));
            directions.push_back(state(3));
        }
        origin(0) = xs.back();
        origin(1) = ys.back();
        origin(2) = yaws.back();
    }

    return {xs, ys, yaws, directions};
}

vector<Path> calc_rs_paths(Vector3d s, Vector3d g, double maxc, double step_size)
{
    vector<Path> paths = generate_path(s, g, maxc, step_size);

    for (Path& path : paths) {
        vector<vector<double>> states = generate_local_course(path.lengths,
                                                         path.ctypes, maxc,
                                                         step_size * maxc);
        for (size_t idx = 0; idx < states[0].size(); ++idx) {
            double ix = states[0][idx];
            double iy = states[1][idx];
            double yaw = states[2][idx];
            int direction = static_cast<int>(states[3][idx]);
            path.x.emplace_back(cos(-s(2)) * ix + sin(-s(2)) * iy + s(0));
            path.y.emplace_back(-sin(-s(2)) * ix + cos(-s(2)) * iy + s(1));
            path.yaw.emplace_back(utils::pi_2_pi(yaw + s(2)));
            path.directions.emplace_back(direction);
        }
    
        for (size_t idx = 0; idx < path.lengths.size(); ++idx) {
            path.lengths[idx] /= maxc;
        }
        path.L /= maxc;
    }

    return paths;
}

Path reeds_shepp_path(Vector3d s, Vector3d g, double maxc, double step_size)
{
    vector<Path> paths = calc_rs_paths(s, g, maxc, step_size);
    int best_path_index = -1;

    for (size_t idx = 0; idx < paths.size(); ++idx) {
        if (best_path_index == -1 ||
            abs(paths[best_path_index].L) > abs(paths[idx].L)) {
            best_path_index = idx;
        }
    }

    if (best_path_index == -1) {
        return Path();
    }

    return paths[best_path_index];
}
