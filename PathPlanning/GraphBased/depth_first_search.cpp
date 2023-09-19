#include <cmath>
#include <stack>
#include <vector>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include <fmt/core.h>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"

using std::stack;
using std::vector;
using std::shared_ptr;
using std::unordered_map;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

class Node
{
public:
    double x;
    double y;
    double cost;
    double parent_index;
    shared_ptr<Node> parent; 

    Node() {}
    Node(double _x, double _y, double _cost, int _parent_index, shared_ptr<Node> _node) {
        x = _x;
        y = _y;
        cost = _cost;
        parent_index = _parent_index;
        parent = _node;
    }
    ~Node() {}
};

class DepthFirstSearchPlanner
{
private:
    double minx;
    double miny;
    double maxx;
    double maxy;
    double xwidth;
    double ywidth;
    double map_resolution;
    double robot_radius;
    vector<vector<bool>> obstacle_map;
    vector<vector<double>> motion;
public:
    DepthFirstSearchPlanner() {}
    DepthFirstSearchPlanner(vector<double> ox, vector<double> oy, double reso, double radius) {
        map_resolution = reso;
        robot_radius = radius;
        calc_obstacle_map(ox, oy);
        motion = get_motion_model();
    }
    ~DepthFirstSearchPlanner() {}
    void calc_obstacle_map(const vector<double>& ox, const vector<double>& oy);
    vector<vector<double>> get_motion_model(void);
    double calc_grid_position(int index, double minp);
    vector<vector<double>> planning(double sx, double sy, double gx, double gy);
    double calc_xyindex(double position, double min_pos);
    double calc_grid_index(shared_ptr<Node> node);
    bool verify_node(shared_ptr<Node> node);
    vector<vector<double>> calc_final_path(shared_ptr<Node> ngoal, 
        unordered_map<double, shared_ptr<Node>>& closed_set);
};

void DepthFirstSearchPlanner::calc_obstacle_map(const vector<double>& ox, const vector<double>& oy)
{
    minx = round(Utils::min(ox));
    miny = round(Utils::min(oy));
    maxx = round(Utils::max(ox));
    maxy = round(Utils::max(oy));
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

vector<vector<double>> DepthFirstSearchPlanner::get_motion_model(void)
{
    vector<vector<double>> motion = {{1, 0, 1}, {0, 1, 1}, {-1, 0, 1}, {0, -1, 1},
                {-1, -1, sqrt(2)}, {-1, 1, sqrt(2)}, {1, -1, sqrt(2)}, {1, 1, sqrt(2)}};
    return motion;
}

double DepthFirstSearchPlanner::calc_grid_position(int index, double minp)
{
    return index * map_resolution + minp;
}


vector<vector<double>> DepthFirstSearchPlanner::planning(double sx, double sy, double gx, double gy)
{
    shared_ptr<Node> nstart = std::make_shared<Node>(calc_xyindex(sx, minx),
                                calc_xyindex(sy, miny), 0.0, -1, nullptr);
    shared_ptr<Node> ngoal = std::make_shared<Node>(calc_xyindex(gx, minx),
                                calc_xyindex(gy, miny), 0.0, -1, nullptr);
    
    unordered_map<double, shared_ptr<Node>> open_set;
    unordered_map<double, shared_ptr<Node>> closed_set;
    open_set[calc_grid_index(nstart)] = nstart;
    stack<shared_ptr<Node>> node_stack;
    node_stack.emplace(nstart);

    while (true) {
        if (open_set.size() == 0 || node_stack.empty()) {
            fmt::print("Open set is empty..\n");
            break;
        }

        shared_ptr<Node> current = node_stack.top();
        node_stack.pop();
        double c_id = calc_grid_index(current);
        open_set.erase(c_id);

        if (show_animation) {
            plt::plot({calc_grid_position(current->x, minx)},
                    {calc_grid_position(current->y, miny)}, "xc");
            plt::pause(0.01);
        }
        if (current->x == ngoal->x && current->y == ngoal->y) {
            fmt::print("Find goal\n");
            ngoal->parent_index = current->parent_index;
            ngoal->cost = current->cost;
            break;
        }

        for (const vector<double>& m : motion) {
            shared_ptr<Node> node = std::make_shared<Node>(current->x + m[0], current->y + m[1],
                                                    current->cost + m[2], c_id, nullptr);
            double n_id = calc_grid_index(node);
            if (!verify_node(node)) {
                continue;
            }

            if (closed_set.find(n_id) == closed_set.end()) {
                node->parent = current;
                open_set[n_id] = node;
                closed_set[n_id] = node;
                node_stack.emplace(node);
            }
        }
    }
    vector<vector<double>> path = calc_final_path(ngoal, closed_set);
    
    return path;
}

double DepthFirstSearchPlanner::calc_xyindex(double position, double min_pos)
{
    return round((position - min_pos) / map_resolution);
}

double DepthFirstSearchPlanner::calc_grid_index(shared_ptr<Node> node)
{
    return (node->y - miny) * xwidth + (node->x - minx);
}

bool DepthFirstSearchPlanner::verify_node(shared_ptr<Node> node)
{
    double px = calc_grid_position(node->x, minx);
    double py = calc_grid_position(node->y, miny);
    if (px < minx || py < miny || px >= maxx || py >= maxy || obstacle_map[node->x][node->y]) {
        return false;
    }

    return true;
}

vector<vector<double>> DepthFirstSearchPlanner::calc_final_path(shared_ptr<Node> ngoal, 
        unordered_map<double, shared_ptr<Node>>& closed_set)
{
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

int main(int argc, char** argv)
{
    double start_x = 10.0;
    double start_y = 10.0;
    double goal_x = 50.0;
    double goal_y = 50.0;
    double grid_size = 2.0;
    double robot_radius = 1.0;

    std::vector<double> obstacle_x;
    std::vector<double> obstacle_y;
    for (int i = -10; i < 60; ++i) {
        obstacle_x.emplace_back(i);
        obstacle_y.emplace_back(-10.0);
    }
    for (int i = -10; i < 60; ++i) {
        obstacle_x.emplace_back(60.0);
        obstacle_y.emplace_back(i);
    }
    for (int i = -10; i < 61; ++i) {
        obstacle_x.emplace_back(i);
        obstacle_y.emplace_back(60.0);
    }
    for (int i = -10; i < 61; ++i) {
        obstacle_x.emplace_back(-10.0);
        obstacle_y.emplace_back(i);
    }
    for (int i = -10; i < 40; ++i) {
        obstacle_x.emplace_back(20.0);
        obstacle_y.emplace_back(i);
    }
    for (int i = 0; i < 40; ++i) {
        obstacle_x.emplace_back(40.0);
        obstacle_y.emplace_back(60.0 - i);
    }
    if (show_animation) {
        plt::plot(obstacle_x, obstacle_y, ".k");
        plt::plot({start_x}, {start_y}, "og");
        plt::plot({goal_x}, {goal_x}, "xb");
        plt::grid(true);
        plt::axis("equal");
    }

    DepthFirstSearchPlanner dfs = DepthFirstSearchPlanner(obstacle_x, obstacle_y, grid_size, robot_radius);
    vector<vector<double>> path = dfs.planning(start_x, start_y, goal_x, goal_y);

    if (show_animation) {
        plt::plot(path[0], path[1], "-r");
        plt::pause(0.01);
        plt::show();
    }

    return 0;
}