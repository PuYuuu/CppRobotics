#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

#include <fmt/core.h>
#include <Eigen/Core>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;
constexpr double KP = 5.0;
constexpr double ETA = 100.0;
constexpr double AREA_WIDTH = 30.0;
constexpr double OSCILLATIONS_DETECTION_LENGTH = 3;

double calc_attractive_potential(Vector2d xy, Vector2d g)
{
    return 0.5 * KP * hypot(xy[0] - g[0], xy[1] - g[1]);
}

double calc_repulsive_potential(Vector2d xy, const vector<vector<double>>& obs, double rr)
{
    if (obs[0].size() < 1) {
        return 0.0;
    }
    int minid = 0;
    double dmin = hypot(xy[0] - obs[0][0], xy[1] - obs[1][0]);
    for (size_t idx = 1; idx < obs.size(); ++idx) {
        double d = hypot(xy[0] - obs[0][idx], xy[1] - obs[1][idx]);
        if (dmin >= d) {
            dmin = d;
            minid = idx;
        }
    }
    double dq = hypot(xy[0] - obs[0][minid], xy[1] - obs[1][minid]);

    if (dq <= rr) {
        if (dq <= 0.1) {
            dq = 0.1;
        }
        return 0.5 * ETA * pow(1.0 / dq - 1.0 / rr, 2); 
    }
        
    return 0.0;
}

Vector2d calc_potential_field(Vector2d s, Vector2d g, const vector<vector<double>>& obs, 
    double reso, double rr, vector<vector<double>>& pmap)
{
    double minx = std::min(Utils::min(obs[0]), std::min(s[0], g[0])) - AREA_WIDTH / 2.0;
    double miny = std::min(Utils::min(obs[1]), std::min(s[1], g[1])) - AREA_WIDTH / 2.0;
    double maxx = std::max(Utils::max(obs[0]), std::max(s[0], g[0])) + AREA_WIDTH / 2.0;
    double maxy = std::max(Utils::max(obs[1]), std::max(s[1], g[1])) + AREA_WIDTH / 2.0;
    int xw = static_cast<int>(round((maxx - minx) / reso));
    int yw = static_cast<int>(round((maxy - miny) / reso));
    vector<vector<double>> pmap_tmp(xw, vector<double>(yw, 0));

    for (int ix = 0; ix < xw; ++ix) {
        double x = ix * reso + minx; 
        for (int iy = 0; iy < yw; ++iy) {
            double y = iy * reso + miny;
            double ug = calc_attractive_potential({x, y}, g);
            double uo = calc_repulsive_potential({x, y}, obs, rr);
            double uf = ug + uo;
            pmap_tmp[ix][iy] = uf;
        }
    }
    pmap = pmap_tmp;

    return {minx, miny};
}

vector<vector<int>> get_motion_model(void)
{
    vector<vector<int>> motion = {{1, 0}, {0, 1}, {-1, 0}, {0, -1},
                                {-1, -1}, {-1, 1}, {1, -1}, {1, 1}};
    return motion;
}

vector<vector<double>> potential_field_planning(Vector2d start, Vector2d goal,
                    const vector<vector<double>>& obs, double reso, double rr)
{
    vector<vector<double>> pmap;
    Vector2d minxy = calc_potential_field(start, goal, obs, reso, rr, pmap);
    double d = hypot(start[0] - goal[0], start[1] - goal[1]);
    int ix = round((start[0] - minxy[0]) / reso);
    int iy = round((start[1] - minxy[1]) / reso);
    int gix = round((goal[0] - minxy[0]) / reso);
    int giy = round((goal[1] - minxy[1]) / reso);

    if (show_animation) {
        int nrows = pmap.size();
        int ncols = pmap[0].size();
        
        vector<float> im(nrows * ncols);
        double max_pmap = -1;
        for (size_t j = 0; j < ncols; ++j) {
            for (size_t i = 0; i < nrows; ++i) {
                im.at(ncols * i + j) = pmap[i][j];
                max_pmap = std::max(max_pmap, pmap[i][j]);
            }
        }

        const float* imptr = &(im[0]);
        const int colors = 1;

        plt::grid(true);
        plt::axis("equal");
        // There is currently no good way to display it, which also makes 
        // debugging difficult, and the method currently fails.
        // plt::imshow(imptr, ncols, nrows, colors);

        plt::plot({static_cast<double>(ix)}, {static_cast<double>(iy)}, "*k");
        plt::plot({static_cast<double>(gix)}, {static_cast<double>(giy)}, "*m");
    }

    vector<vector<double>> path = {{start[0]}, {start[1]}};
    vector<vector<int>> motions = get_motion_model();
    
    while (d >= reso) {
        double minp = std::numeric_limits<double>::max() - 10;
        int minix = -1;
        int miniy = -1;
        for (vector<int> motion : motions) {
            int inx = ix + motion[0];
            int iny = iy + motion[1];
            double p = std::numeric_limits<double>::max() - 10;
            if (inx < pmap.size() && iny < pmap[0].size() && inx >= 0 && iny >= 0) {
                p = pmap[inx][iny];
            }

            if (minp > p) {
                minp = p;
                minix = inx;
                miniy = iny;
            }
        }
        ix = minix;
        iy = miniy;
        double xp = ix * reso + minxy[0];
        double yp = iy * reso + minxy[1];
        d = hypot(goal[0] - xp, goal[1] - yp);
        path[0].push_back(xp);
        path[1].push_back(yp);

        if (show_animation) {
            plt::plot({static_cast<double>(ix)}, {static_cast<double>(iy)}, ".r");
            plt::pause(0.01);
        }
    }

    return path;
}

int main (int argc, char** argv)
{
    Vector2d start(0, 10);
    Vector2d goal(30, 30);
    double grid_size = 0.5;
    double robot_radius = 5.0;
    vector<vector<double>> obstacles = {{15.0, 5.0, 20.0, 25.0}, 
                                        {25.0, 15.0, 26.0, 25.0}};

    potential_field_planning(start, goal, obstacles, grid_size, robot_radius);
    if (show_animation) {
        plt::show();
    }

    return 0;
}
