#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

#include <fmt/core.h>
#include <Eigen/Core>

#include "utils.hpp"
#include "cruise_line.hpp"
#include "cubic_spline.hpp"
#include "matplotlibcpp.h"
#include "quintic_polynomial.hpp"
#include "quartic_polynomial.hpp"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr double ROAD_WIDTH = 8.0;
constexpr double ROAD_SAMPLE_STEP = 1.0;
constexpr double TARGET_SPEED = 30.0 / 3.6;
constexpr double SPEED_SAMPLE_STEP = 5.0 / 3.6;

constexpr double T_STEP = 0.15;
constexpr double K_JERK = 0.1;
constexpr double K_TIME = 1.0;
constexpr double K_V_DIFF = 1.0;
constexpr double K_OFFSET = 1.5;
constexpr double K_COLLISION = 500;

constexpr double MAX_SPEED = 50.0 / 3.6;
constexpr double MAX_ACCEL = 8.0;
constexpr double MAX_CURVATURE = 6.0;

class Path
{
public:
    vector<double> t;
    double cost = 0.0;

    vector<double> l;
    vector<double> l_v;
    vector<double> l_a;
    vector<double> l_jerk;

    vector<double> s;
    vector<double> s_v;
    vector<double> s_a;
    vector<double> s_jerk;

    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    vector<double> ds;
    vector<double> curv;

    Path() {}
    ~Path() {}

    void SL_2_XY(CubicSpline2D& ref_path);
    void calc_yaw_curv(void);
    bool operator<(const Path& other) const {
        return cost < other.cost;
    }
};

void Path::SL_2_XY(CubicSpline2D& ref_path)
{
    x.clear();
    y.clear();

    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] > ref_path.s.back()) {
            break;
        }

        Vector2d xy_ref = ref_path.calc_position(s[i]);
        double yaw =ref_path.calc_yaw(s[i]);
        double x_ref = xy_ref[0] + l[i] * cos(yaw + M_PI_2);
        double y_ref = xy_ref[1] + l[i] * sin(yaw + M_PI_2);

        x.push_back(x_ref);
        y.push_back(y_ref);
    }
}

void Path::calc_yaw_curv(void)
{
    yaw.clear();
    curv.clear();
    ds.clear();
    
    for (size_t i = 0; i + 1 < x.size(); ++i) {
        double dx = x[i + 1] - x[i];
        double dy = y[i + 1] - y[i];
        ds.push_back(hypot(dx, dy));
        yaw.push_back(atan2(dy, dx));
    }

    if (yaw.empty()) {
        return ;
    }
    yaw.push_back(yaw.back());
    ds.push_back(ds.back());

    for (size_t i = 0; i + 1 < yaw.size(); ++i) {
        curv.emplace_back((yaw[i + 1] - yaw[i]) / ds[i]);
    }
}

vector<vector<double>> get_reference_line(vector<double> cx, vector<double> cy, CubicSpline2D& spline)
{
    vector<double> x;
    vector<double> y;
    for (size_t idx = 0; idx < cx.size(); idx += 3) {
        x.push_back(cx[idx]);
        y.push_back(cy[idx]);
    }

    vector<vector<double>> traj = CubicSpline2D::calc_spline_course(x, y, 0.1);
    spline = CubicSpline2D(x, y);
    
    return traj;
}

bool verify_path(const Path& path)
{
    for (size_t i = 0; i < path.s_v.size(); ++i) {
        if (path.s_v[i] > MAX_SPEED || abs(path.s_a[i]) > MAX_ACCEL ||
            abs(path.curv[i]) > MAX_CURVATURE) {
            return false;
        }
    }

    return true;
}

double is_path_collision(const Path& path, const utils::VehicleConfig& vc, const vector<vector<double>>& obs)
{
    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    for (size_t i = 0; i < path.x.size(); i += 3) {
        x.push_back(path.x[i]);
        y.push_back(path.y[i]);
        yaw.push_back(path.yaw[i]);
    }

    for (size_t i = 0; i < x.size(); ++i) {
        double d = 1.8;
        double dl = (vc.RF - vc.RB) / 2.0;
        double r = hypot((vc.RF + vc.RB) / 2.0, vc.W / 2.0) + d;

        double cx = x[i] + dl * cos(yaw[i]);
        double cy = y[i] + dl * sin(yaw[i]);

        for (size_t j = 0; j < obs[0].size(); ++j) {
            double xo = obs[0][j] - cx;
            double yo = obs[1][j] - cy;
            double dx = xo * cos(yaw[i]) + yo * sin(yaw[i]);
            double dy = -xo * sin(yaw[i]) + yo * cos(yaw[i]);

            if (abs(dx) < r && abs(dy) < vc.W / 2 + d) {
                return 1.0;
            }
        }
    }

    return 0.0;
}

vector<Path> sampling_paths(
    double l0, double l0_v, double l0_a, double s0, double s0_v, double s0_a,
    CubicSpline2D& ref_path, const utils::VehicleConfig& vc, const vector<vector<double>>& obs)
{
    vector<Path> paths;

    for (double s1_v = TARGET_SPEED * 0.6; s1_v < TARGET_SPEED * 1.4; s1_v += TARGET_SPEED * 0.2) {
        for (double t1 = 4.5; t1 < 5.5; t1 += 0.2) {
            Path path_pre;
            QuarticPolynomial path_lon(s0, s0_v, s0_a, s1_v, 0.0, t1);
            
            for (double t = 0.0; t < t1; t += T_STEP) {
                path_pre.t.push_back(t);
                path_pre.s.push_back(path_lon.calc_point(t));
                path_pre.s_v.push_back(path_lon.calc_first_derivative(t));
                path_pre.s_a.push_back(path_lon.calc_second_derivative(t));
                path_pre.s_jerk.push_back(path_lon.calc_third_derivative(t));
            }
            
            for (double l1 = -ROAD_WIDTH; l1 < ROAD_WIDTH; l1 += ROAD_SAMPLE_STEP) {
                Path path = path_pre;
                QuinticPolynomial path_lat(l0, l0_v, l0_a, l1, 0.0, 0.0, t1);

                for (double t : path_pre.t) {
                    path.l.push_back(path_lat.calc_point(t));
                    path.l_v.push_back(path_lat.calc_first_derivative(t));
                    path.l_a.push_back(path_lat.calc_second_derivative(t));
                    path.l_jerk.push_back(path_lat.calc_third_derivative(t));
                }

                path.SL_2_XY(ref_path);
                path.calc_yaw_curv();
                if (path.yaw.empty()) {
                    continue;
                }

                double l_jerk_sum = 0.0;
                double s_jerk_sum = 0.0;
                double v_diff = abs(TARGET_SPEED - path.s_v.back());
                for (size_t i = 0; i < path.l_jerk.size(); ++i) {
                    l_jerk_sum += abs(path.l_jerk[i]);
                    s_jerk_sum += abs(path.s_jerk[i]);
                }
                
                path.cost = K_JERK * (l_jerk_sum + s_jerk_sum) +
                        K_V_DIFF * v_diff + K_TIME * t1 * 2 +
                        K_OFFSET * abs(path.l.back()) +
                        K_COLLISION * is_path_collision(path, vc, obs);

                paths.emplace_back(path);
            }
        }
    }

    return paths;
}

Path extract_optimal_path(vector<Path>& paths)
{
    Path path;
    if (paths.empty()) {
        return path;
    }

    std::sort(paths.begin(), paths.end());
    for (Path& p : paths) {
        if (verify_path(p)) {
            path = p;
            return path;
        }
    }

    return paths.back();
}

Path lattice_planner(
    double l0, double l0_v, double l0_a, double s0, double s0_v, double s0_a,
    CubicSpline2D& ref_path, const utils::VehicleConfig& vc, const vector<vector<double>>& obs)
{
    vector<Path> paths = sampling_paths(l0, l0_v, l0_a, s0, s0_v, s0_a, ref_path, vc, obs);
    Path path = extract_optimal_path(paths);

    return path;
}

int main()
{
    CruiseLine cruise_line;
    vector<vector<double>> wxy = cruise_line.design_reference_line();
    vector<vector<double>> inxy = cruise_line.design_boundary_in();
    vector<vector<double>> outxy = cruise_line.design_boundary_out();
    vector<vector<double>> obs = {{50, 96, 70, 40, 25},{10, 25, 40, 50, 75}};
    CubicSpline2D spline;
    vector<vector<double>> traj = get_reference_line(wxy[0], wxy[1], spline);

    utils::VehicleConfig vc;
    vc.RF = 6.75;
    vc.RB = 1.5;
    vc.W = 4.5;
    vc.WD = 0.7 * vc.W;
    vc.WB = 5.25;
    vc.TR = 0.75;
    vc.TW = 1.5;

    double l0 = 0.0;    // current lateral position [m]
    double l0_v = 0.0;  // current lateral speed [m/s]
    double l0_a = 0.0;  // current lateral acceleration [m/s]
    double s0 = 0.0;    // current course position
    double s0_v = 30.0 / 3.6;   // current speed [m/s]
    double s0_a = 0.0;

    while (true) {
        Path path = lattice_planner(l0, l0_v, l0_a, s0, s0_v, s0_a, spline, vc, obs);

        if (path.x.empty()) {
            fmt::print("No feasible path found!!\n");
            break;
        }

        l0 = path.l[1];
        l0_v = path.l_v[1];
        l0_a = path.l_a[1];
        s0 = path.s[1];
        s0_v = path.s_v[1];
        s0_a = path.s_a[1];

        if (hypot(path.x[1] - traj[0].back(), path.y[1] - traj[1].back()) <= 2.0) {
            fmt::print("Goal\n");
            break;
        }

        double dy = (path.yaw[2] - path.yaw[1]) / path.ds[1];
        double steer = utils::pi_2_pi(atan(1.2 * vc.WB * dy));
        plt::cla();
        plt::plot(wxy[0], wxy[1], {{"linestyle", "--"}, {"color", "gray"}});
        plt::plot(inxy[0], inxy[1], {{"linewidth", "2"}, {"color", "k"}});
        plt::plot(outxy[0], outxy[1], {{"linewidth", "2"}, {"color", "k"}});
        plt::plot(path.x, path.y, "-r");
        plt::plot(obs[0], obs[1], "ok");
        utils::draw_vehicle({path.x[1], path.y[1], path.yaw[1]}, steer, vc);
        plt::title("Lattice Planner in Cruising Scene V[km/h]:" +
                        std::to_string(s0_v * 3.6).substr(0,4));
        plt::axis("equal");
        plt::pause(0.0001);
    }
    plt::show();

    return 0;
}
