#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include <fmt/core.h>

#include "utils.hpp"
#include "matplotlibcpp.h"
#include "GraphSearchPlanner.hpp"

using std::vector;
using std::shared_ptr;
using std::unordered_map;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

class BidirectionalAStarPlanner : public GraphSearchPlanner
{
private:
    shared_ptr<Node> get_mincost_node(
        const unordered_map<double, shared_ptr<Node>>& node_set, shared_ptr<Node> goal);
    double calc_heuristic(shared_ptr<Node> node1, shared_ptr<Node> node2);
    vector<vector<double>> calc_final_bidirectional_path(shared_ptr<Node> n1, shared_ptr<Node> n2,
        unordered_map<double, shared_ptr<Node>>& set1, unordered_map<double, shared_ptr<Node>>& set2);
public:
    BidirectionalAStarPlanner() {}
    BidirectionalAStarPlanner(vector<double> ox, vector<double> oy, double reso, double radius) :
        GraphSearchPlanner(ox, oy, reso, radius) {}
    ~BidirectionalAStarPlanner() override {}

    vector<vector<double>> planning(double sx, double sy, double gx, double gy) override;
};

vector<vector<double>> BidirectionalAStarPlanner::planning(double sx, double sy, double gx, double gy)
{
    shared_ptr<Node> nstart = std::make_shared<Node>(calc_xyindex(sx, get_minx()),
                                calc_xyindex(sy, get_miny()), 0.0, -1, nullptr);
    shared_ptr<Node> ngoal = std::make_shared<Node>(calc_xyindex(gx, get_minx()),
                                calc_xyindex(gy, get_miny()), 0.0, -1, nullptr);
    unordered_map<double, shared_ptr<Node>> open_set_A, closed_set_A;
    unordered_map<double, shared_ptr<Node>> open_set_B, closed_set_B;
    open_set_A[calc_grid_index(nstart)] = nstart;
    open_set_B[calc_grid_index(ngoal)] = ngoal;
    shared_ptr<Node> current_A = nstart;
    shared_ptr<Node> current_B = ngoal;
    shared_ptr<Node> meet_point_A = nullptr;
    shared_ptr<Node> meet_point_B = nullptr;

    while (true) {
        if (open_set_A.size() == 0 || open_set_B.size() == 0) {
            fmt::print("Open set is empty..\n");
            break;
        }

        current_A = get_mincost_node(open_set_A, current_B);
        double c_id_A = calc_grid_index(current_A);
        current_B = get_mincost_node(open_set_B, current_A);
        double c_id_B = calc_grid_index(current_B);
        open_set_A.erase(c_id_A);
        open_set_B.erase(c_id_B);

        if (show_animation) {
            plt::plot({calc_grid_position(current_A->x, get_minx())},
                    {calc_grid_position(current_A->y, get_miny())}, "xc");
            plt::plot({calc_grid_position(current_B->x, get_minx())},
                    {calc_grid_position(current_B->y, get_miny())}, "xy");
            if (closed_set_A.size() % 10 == 0) {
                plt::pause(0.001);
            }
        }
        if (current_A->x == current_B->x && current_A->y == current_B->y) {
            fmt::print("Find goal\n");
            meet_point_A = current_A;
            meet_point_B = current_B;
            break;
        }

        closed_set_A[c_id_A] = current_A;
        closed_set_B[c_id_B] = current_B;
        for (const vector<double>& m : get_motion()) {
            vector<shared_ptr<Node>> nodes = {  std::make_shared<Node>(current_A->x + m[0], current_A->y + m[1],
                                                    current_A->cost + m[2], c_id_A, current_A),
                                                std::make_shared<Node>(current_B->x + m[0], current_B->y + m[1],
                                                    current_B->cost + m[2], c_id_B, current_B)};
            vector<double> n_ids = {calc_grid_index(nodes[0]), calc_grid_index(nodes[1])};

            bool continue_A = (!verify_node(nodes[0]) || closed_set_A.find(n_ids[0]) != closed_set_A.end());
            bool continue_B = (!verify_node(nodes[1]) || closed_set_B.find(n_ids[1]) != closed_set_B.end());
            
            if (!continue_A && 
                (open_set_A.find(n_ids[0]) == open_set_A.end() || open_set_A[n_ids[0]]->cost >= nodes[0]->cost)) {
                open_set_A[n_ids[0]] = nodes[0];
            }
            if (!continue_B && 
                (open_set_B.find(n_ids[1]) == open_set_B.end() || open_set_B[n_ids[1]]->cost >= nodes[1]->cost)) {
                open_set_B[n_ids[1]] = nodes[1];
            }
        }
    }
    vector<vector<double>> path = calc_final_bidirectional_path(
                        meet_point_A, meet_point_B, closed_set_A, closed_set_B);

    return path;
}

shared_ptr<Node> BidirectionalAStarPlanner::get_mincost_node(
    const unordered_map<double, shared_ptr<Node>>& node_set, shared_ptr<Node> goal)
{
    shared_ptr<Node> min_node = nullptr;
    for (const auto& n : node_set ) {
        if (min_node == nullptr || 
            ((min_node->cost + calc_heuristic(min_node, goal)) > 
            (n.second->cost + calc_heuristic(n.second, goal)))) {
            min_node = n.second;
        }
    }
    return min_node;
}

double BidirectionalAStarPlanner::calc_heuristic(shared_ptr<Node> node1, shared_ptr<Node> node2)
{
    double weight = 1.0;
    double distance = weight * hypot(node1->x - node2->x, node1->y - node2->y);

    return distance;
}

vector<vector<double>> BidirectionalAStarPlanner::calc_final_bidirectional_path(
    shared_ptr<Node> n1, shared_ptr<Node> n2,
    unordered_map<double, shared_ptr<Node>>& set1, unordered_map<double, shared_ptr<Node>>& set2)
{
    vector<vector<double>> path1 = calc_final_path(n1, set1);
    vector<vector<double>> path2 = calc_final_path(n2, set2);
    std::reverse(path1[0].begin(), path1[0].end());
    std::reverse(path1[1].begin(), path1[1].end());
    path1[0].insert(path1[0].end(), path2[0].begin(), path2[0].end());
    path1[1].insert(path1[1].end(), path2[1].begin(), path2[1].end());

    return path1;
}

int main(int argc, char** argv)
{
    double start_x = 10;
    double start_y = 10;
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
        plt::plot({goal_x}, {goal_x}, "ob");
        plt::title("Bidirectional-A*");
        plt::grid(true);
        plt::axis("equal");
    }

    BidirectionalAStarPlanner bidir_a_star(obstacle_x, obstacle_y, grid_size, robot_radius);
    vector<vector<double>> path = bidir_a_star.planning(start_x, start_y, goal_x, goal_y);

    if (show_animation) {
        plt::plot(path[0], path[1], "-r");
        plt::pause(0.01);
        plt::show();
    }

    return 0;
}
