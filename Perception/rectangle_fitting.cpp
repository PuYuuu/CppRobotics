#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <tuple>
#include <set>

#include <fmt/core.h>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include "simulator.hpp"
#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::set;
using std::vector;
using std::tuple;
using namespace Eigen;
namespace plt = matplotlibcpp;
constexpr double DT = 0.2;
constexpr double SIM_TIME = 30.0;
constexpr bool show_animation = true;

enum Criteria {AREA, CLOSENESS, VARIANCE};

class RectangleData
{
public:
    vector<double> a;
    vector<double> b;
    vector<double> c;
    vector<double> rect_c_x;
    vector<double> rect_c_y;

    RectangleData() {
        a.resize(4);
        b.resize(4);
        c.resize(4);
        rect_c_x.resize(5);
        rect_c_y.resize(5);
    }
    ~RectangleData() {}

    void plot(void);
    void calc_rect_contour(void);
    tuple<double, double> calc_cross_point(
        vector<double> _a, vector<double> _b, vector<double> _c);
};

void RectangleData::plot(void)
{
    calc_rect_contour();
    plt::plot(rect_c_x, rect_c_y, "-r");
}

void RectangleData::calc_rect_contour(void)
{
    std::tie(rect_c_x[0], rect_c_y[0]) = calc_cross_point(
        {a[0], a[1]}, {b[0], b[1]}, {c[0], c[1]});
    std::tie(rect_c_x[1], rect_c_y[1]) = calc_cross_point(
        {a[1], a[2]}, {b[1], b[2]}, {c[1], c[2]});
    std::tie(rect_c_x[2], rect_c_y[2]) = calc_cross_point(
        {a[2], a[3]}, {b[2], b[3]}, {c[2], c[3]});
    std::tie(rect_c_x[3], rect_c_y[3]) = calc_cross_point(
        {a[3], a[0]}, {b[3], b[0]}, {c[3], c[0]});
    rect_c_x[4] = rect_c_x[0];
    rect_c_y[4] = rect_c_y[0];
}

tuple<double, double> RectangleData::calc_cross_point(
        vector<double> _a, vector<double> _b, vector<double> _c)
{
    double x = (_b[0] * -_c[1] - _b[1] * - _c[0]) / (_a[0] * _b[1] - _a[1] * _b[0]);
    double y = (_a[1] * -_c[0] - _a[0] * - _c[1]) / (_a[0] * _b[1] - _a[1] * _b[0]);
    
    return std::make_tuple(x, y);
}

class LShapeFitting
{
private:
    Criteria criteria;
    double min_dist_of_closeness_criteria;
    double d_theta_deg_for_search;
    double R0 = 3.0;
    double Rd = 0.001;

    Vector4d find_min_max(const vector<double>& c1, const vector<double>& c2);
    vector<vector<int>> adoptive_range_segmentation(const vector<vector<double>>& oxy);
    RectangleData rectangle_search(const vector<vector<double>>& cxy);
    double calc_area_criterion(const vector<double>& c1, const vector<double>& c2);
    double calc_closeness_criterion(const vector<double>& c1, const vector<double>& c2);
    double calc_variance_criterion(const vector<double>& c1, const vector<double>& c2);
public:
    LShapeFitting() {
        criteria = VARIANCE;
        min_dist_of_closeness_criteria = 0.01;
        d_theta_deg_for_search = 1.0;
        R0 = 3.0;
        Rd = 0.001;
    }
    ~LShapeFitting() {}

    vector<RectangleData> fitting(const vector<vector<double>>& oxy, 
                                        vector<vector<int>>& id_sets);
};

Vector4d LShapeFitting::find_min_max(const vector<double>& c1, const vector<double>& c2)
{
    double c1_max = utils::max(c1);
    double c2_max = utils::max(c2);
    double c1_min = utils::min(c1);
    double c2_min = utils::min(c2);
    Vector4d min_max;
    min_max << c1_max, c1_min, c2_max, c2_min;

    return min_max;
}

vector<vector<int>> LShapeFitting::adoptive_range_segmentation(
                                    const vector<vector<double>>& oxy)
{
    vector<set<int>> segment_list;

    for (size_t i = 0; i < oxy[0].size(); ++i) {
        set<int> c;
        double r = R0 + Rd * sqrt(oxy[0][i] * oxy[0][i] + oxy[1][i] * oxy[1][i]);
        for (size_t j = 0; j < oxy[0].size(); ++j) {
            double d = hypot(oxy[0][i] - oxy[0][j], oxy[1][i] - oxy[1][j]);
            if (d <= r) {
                c.insert(j);
            }
        }
        segment_list.emplace_back(c);
    }

    for (size_t i = 0; i < segment_list.size() - 1; ++i) {
        for (size_t j = i + 1; j < segment_list.size(); ++j) {
            std::set<int> tmp;
            std::set_intersection(segment_list[i].begin(), segment_list[i].end(),
                                    segment_list[j].begin(), segment_list[j].end(), 
                                    std::inserter(tmp, tmp.begin()));
            if (!tmp.empty()) {
                tmp.clear();
                std::set_union(segment_list[i].begin(), segment_list[i].end(),
                                segment_list[j].begin(), segment_list[j].end(), 
                                std::inserter(tmp, tmp.begin()));
                segment_list[i].clear();
                segment_list[j].clear();
                segment_list[i] = tmp;
            }
        }
    }
    
    vector<vector<int>> id_sets;
    for (size_t idx = 0; idx < segment_list.size(); ++idx) {
        if (!segment_list[idx].empty()) {
            vector<int> tmp;
            for (int id : segment_list[idx]) {
                tmp.push_back(id);
            }
            id_sets.emplace_back(tmp);
        }
    }

    return id_sets;
}

double LShapeFitting::calc_area_criterion(const vector<double>& c1, const vector<double>& c2)
{
    Vector4d min_max = find_min_max(c1, c2);
    double alpha = -(min_max(0) - min_max(1)) * (min_max(0) - min_max(1));
    
    return alpha;
}

double LShapeFitting::calc_closeness_criterion(const vector<double>& c1, const vector<double>& c2)
{
    Vector4d min_max = find_min_max(c1, c2);
    double beta = 0.0;
    for (size_t idx = 0; idx < c1.size(); ++idx) {
        double d1 = std::min(min_max(0) - c1[idx], c1[idx] - min_max(1));
        double d2 = std::min(min_max(2) - c2[idx], c2[idx] - min_max(3));
        double d = std::min(std::min(d1, d2), min_dist_of_closeness_criteria);
        beta += (1.0 / d);
    }

    return beta;
}

double LShapeFitting::calc_variance_criterion(const vector<double>& c1, const vector<double>& c2)
{
    Vector4d min_max = find_min_max(c1, c2);
    vector<double> e1;
    vector<double> e2;
    for (size_t idx = 0; idx < c1.size(); ++idx) {
        double d1 = std::min(min_max(0) - c1[idx], c1[idx] - min_max(1));
        double d2 = std::min(min_max(2) - c2[idx], c2[idx] - min_max(3));
        if (d1 < d2) {
            e1.emplace_back(d1);
        } else {
            e2.emplace_back(d2);
        }
    }

    double v1 = -utils::variance(e1);
    double v2 = -utils::variance(e2);
    double gamma = v1 + v2;

    return gamma;
}

RectangleData LShapeFitting::rectangle_search(const vector<vector<double>>& cxy)
{
    double d_theta = d_theta_deg_for_search * M_PI / 180.0;
    Vector2d min_cost(-1e6, 0.);

    for (double theta = 0; theta < M_PI_2 - d_theta; theta += d_theta) {
        vector<double> c1;
        vector<double> c2;
        Matrix2d l_rot = utils::rotation_matrix2d(theta);
        for (size_t idx = 0; idx < cxy[0].size(); ++idx) {
            Vector2d xy(cxy[0][idx], cxy[1][idx]);
            Matrix<double, 1, 2> le_xy = xy.transpose() * l_rot;
            c1.push_back(le_xy(0, 0));
            c2.push_back(le_xy(0, 1));
        }
        double cost = 0.0;
        if (criteria == AREA) {
            cost = calc_area_criterion(c1, c2);
        } else if (criteria == CLOSENESS) {
            cost = calc_closeness_criterion(c1, c2);
        } else if (criteria == VARIANCE) {
            cost = calc_variance_criterion(c1, c2);
        }

        if (min_cost(0) < cost) {
            min_cost(0) = cost;
            min_cost(1) = theta;
        }
    }

    double sin_s = sin(min_cost(1));
    double cos_s = cos(min_cost(1));
    vector<double> c1_s;
    vector<double> c2_s;

    for (size_t idx = 0; idx < cxy[0].size(); ++idx) {
        c1_s.emplace_back(cxy[0][idx] * cos_s + cxy[1][idx] * sin_s);
        c2_s.emplace_back(-cxy[0][idx] * sin_s + cxy[1][idx] * cos_s);
    }

    RectangleData rect;
    rect.a[0] = cos_s;
    rect.b[0] = sin_s;
    rect.c[0] = utils::min(c1_s);
    rect.a[1] = -sin_s;
    rect.b[1] = cos_s;
    rect.c[1] = utils::min(c2_s);
    rect.a[2] = cos_s;
    rect.b[2] = sin_s;
    rect.c[2] = utils::max(c1_s);
    rect.a[3] = -sin_s;
    rect.b[3] = cos_s;
    rect.c[3] = utils::max(c2_s);

    return rect;
}

vector<RectangleData> LShapeFitting::fitting(
    const vector<vector<double>>& oxy, vector<vector<int>>& id_sets)
{
    id_sets = adoptive_range_segmentation(oxy);
    vector<RectangleData> rects;

    for (const vector<int>& ids : id_sets) {
        vector<vector<double>> cxy(2);
        for (int id : ids) {
            cxy[0].push_back(oxy[0][id]);
            cxy[1].push_back(oxy[1][id]);
        }
        RectangleData rect = rectangle_search(cxy);
        rects.emplace_back(rect);
    }

    return rects;
}

int main(int argc, char** argv)
{
    double angle_resolution = 3 * M_PI / 180.0;
    VehicleSimulator v1(-10.0, 0.0, M_PI_2, 0.0, 50.0 / 3.6, 3.0, 5.0);
    VehicleSimulator v2(20.0, 10.0, M_PI, 0.0, 50.0 / 3.6, 4.0, 10.0);

    LShapeFitting l_shape_fitting;
    LidarSimulator lidar_sim;
    double time = 0.0;

    while (time <= SIM_TIME) {
        time += DT;

        v1.update(DT, 0.1, 0.0);
        v2.update(DT, 0.1, -0.05);

        vector<vector<int>> id_sets;
        vector<vector<double>> oxy =
            lidar_sim.get_observation_points({v1, v2}, angle_resolution);
        vector<RectangleData> rects = l_shape_fitting.fitting(oxy, id_sets);
        if (show_animation) {
            plt::cla();
            plt::axis("equal");
            plt::plot({0.0}, {0.0}, "*r");
            v1.plot();
            v2.plot();

            for (const vector<int>& ids : id_sets) {
                vector<vector<double>> sets_xy(2);
                for (int id : ids) {
                    sets_xy[0].push_back(oxy[0][id]);
                    sets_xy[1].push_back(oxy[1][id]);
                }
                for (size_t idx = 0; idx < sets_xy[0].size(); ++idx) {
                    plt::plot({0.0, sets_xy[0][idx]}, {0.0, sets_xy[1][idx]}, "-og");
                }
                plt::plot(sets_xy[0], sets_xy[1], "o");
            }

            for (RectangleData rect : rects) {
                rect.plot();
            }

            plt::title("Rectangle Fitting");
            plt::legend({{"loc", "upper right"}});
            plt::pause(0.1);
        }
    }

    return 0;
}
