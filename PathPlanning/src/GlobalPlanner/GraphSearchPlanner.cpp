#include "GraphSearchPlanner.hpp"

using std::shared_ptr;
using std::unordered_map;
using std::vector;

void GraphSearchPlanner::calc_obstacle_map(const vector<double>& ox, const vector<double>& oy) {
    minx = round(utils::min(ox));
    miny = round(utils::min(oy));
    maxx = round(utils::max(ox));
    maxy = round(utils::max(oy));
    fmt::print("min_x: {}\n", minx);
    fmt::print("min_y: {}\n", miny);
    fmt::print("max_x: {}\n", maxx);
    fmt::print("max_y: {}\n", maxy);

    xwidth = round((maxx - minx) / map_resolution);
    ywidth = round((maxy - miny) / map_resolution);
    fmt::print("x_width: {}\n", xwidth);
    fmt::print("y_width: {}\n", ywidth);

    vector<vector<bool>> obsmap(xwidth, vector<bool>(ywidth, false));
    obstacle_map = obsmap;
    for (int ix = 0; ix < xwidth; ++ix) {
        double x = calc_grid_position(ix, minx);
        for (int iy = 0; iy < ywidth; ++iy) {
            double y = calc_grid_position(iy, miny);
            for (int idx = 0; idx < ox.size(); ++idx) {
                double rho = hypot(ox[idx] - x, oy[idx] - y);
                if (rho <= robot_radius) {
                    obstacle_map[ix][iy] = true;
                    break;
                }
            }
        }
    }
}

vector<vector<double>> GraphSearchPlanner::get_motion_model(void) {
    vector<vector<double>> motion = {{1, 0, 1},        {0, 1, 1},         {-1, 0, 1},
                                     {0, -1, 1},       {-1, -1, sqrt(2)}, {-1, 1, sqrt(2)},
                                     {1, -1, sqrt(2)}, {1, 1, sqrt(2)}};
    return motion;
}

double GraphSearchPlanner::calc_grid_position(int index, double minp) {
    return index * map_resolution + minp;
}

double GraphSearchPlanner::calc_xyindex(double position, double min_pos) {
    return round((position - min_pos) / map_resolution);
}

double GraphSearchPlanner::calc_grid_index(shared_ptr<Node> node) {
    return (node->y - miny) * xwidth + (node->x - minx);
}

bool GraphSearchPlanner::verify_node(shared_ptr<Node> node) {
    double px = calc_grid_position(node->x, minx);
    double py = calc_grid_position(node->y, miny);
    if (px < minx || py < miny || px >= maxx || py >= maxy || obstacle_map[node->x][node->y]) {
        return false;
    }

    return true;
}

vector<vector<double>> GraphSearchPlanner::calc_final_path(
    shared_ptr<Node> ngoal, unordered_map<double, shared_ptr<Node>>& closed_set) {
    vector<vector<double>> path = {{calc_grid_position(ngoal->x, minx)},
                                   {calc_grid_position(ngoal->y, miny)}};
    shared_ptr<Node> node = closed_set[ngoal->parent_index];
    while (node != nullptr) {
        path[0].emplace_back(calc_grid_position(node->x, minx));
        path[1].emplace_back(calc_grid_position(node->y, miny));
        node = node->parent;
    }

    return path;
}
