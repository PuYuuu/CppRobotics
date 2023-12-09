#include <cmath>
#include <stack>
#include <vector>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include <fmt/core.h>

#include "utils.hpp"
#include "matplotlibcpp.h"
#include "GraphSearchPlanner.hpp"

using std::stack;
using std::vector;
using std::shared_ptr;
using std::unordered_map;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

class DepthFirstSearchPlanner : public GraphSearchPlanner
{
public:
    DepthFirstSearchPlanner() {}
    DepthFirstSearchPlanner(vector<double> ox, vector<double> oy, double reso, double radius) :
        GraphSearchPlanner(ox, oy, reso, radius) {}
    ~DepthFirstSearchPlanner() {}

    vector<vector<double>> planning(double sx, double sy, double gx, double gy) override;
};


vector<vector<double>> DepthFirstSearchPlanner::planning(double sx, double sy, double gx, double gy)
{
    shared_ptr<Node> nstart = std::make_shared<Node>(calc_xyindex(sx, get_minx()),
                                calc_xyindex(sy, get_miny()), 0.0, -1, nullptr);
    shared_ptr<Node> ngoal = std::make_shared<Node>(calc_xyindex(gx, get_minx()),
                                calc_xyindex(gy, get_miny()), 0.0, -1, nullptr);
    
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
            plt::plot({calc_grid_position(current->x, get_minx())},
                    {calc_grid_position(current->y, get_miny())}, "xc");
            plt::pause(0.01);
        }
        if (current->x == ngoal->x && current->y == ngoal->y) {
            fmt::print("Find goal\n");
            ngoal->parent_index = current->parent_index;
            ngoal->cost = current->cost;
            break;
        }

        for (const vector<double>& m : get_motion()) {
            shared_ptr<Node> node = std::make_shared<Node>(current->x + m[0], current->y + m[1],
                                                    current->cost + m[2], c_id, current);
            double n_id = calc_grid_index(node);
            if (!verify_node(node)) {
                continue;
            }

            if (closed_set.find(n_id) == closed_set.end()) {
                // node->parent = current;
                open_set[n_id] = node;
                closed_set[n_id] = node;
                node_stack.emplace(node);
            }
        }
    }
    vector<vector<double>> path = calc_final_path(ngoal, closed_set);

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
        plt::plot(obstacle_x, obstacle_y, "sk");
        plt::plot({start_x}, {start_y}, "og");
        plt::plot({goal_x}, {goal_x}, "xb");
        plt::grid(true);
        plt::title("Depth-first Search");
        plt::axis("equal");
    }

    DepthFirstSearchPlanner dfs(obstacle_x, obstacle_y, grid_size, robot_radius);
    vector<vector<double>> path = dfs.planning(start_x, start_y, goal_x, goal_y);

    if (show_animation) {
        plt::plot(path[0], path[1], "-r");
        plt::pause(0.01);
        plt::show();
    }

    return 0;
}
