#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <unordered_map>

#include <fmt/core.h>
#include <Eigen/Core>

#include "utils.hpp"
#include "matplotlibcpp.h"

using std::vector;
using std::string;
using namespace Eigen;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

class Node
{
public:
    double x;
    double y;
    double cost;
    int parent_index;

    Node() {}
    Node(double _x, double _y, double _cost = 0., int parent = -1) :
        x(_x), y(_y), cost(_cost), parent_index(parent) {}
    ~Node() {}
};

class PRM
{
private:
    Vector2d start;
    Vector2d goal;
    double robot_radius;
    double min_rand;
    double max_rand;
    int sample_num;
    vector<vector<double>> obstacle_list;
    vector<vector<double>> boundary;
    std::mt19937 engine;

    bool check_collision(Vector2d node);
    bool check_collision(Vector2d node1, Vector2d node2);
    vector<vector<double>> calc_final_path( const vector<vector<double>>& vertex,
                                            const vector<vector<double>>& edge);
    void plot_road_map( const vector<vector<double>>& vertex,
                        const vector<vector<double>>& edge);
    void plot_circle(double x, double y, double size,
                    bool is_fill = false, string style = "-b");
public:
    PRM(Vector2d _start, Vector2d _goal, vector<vector<double>> obs,
        Vector2d rand_area, double _robot_radius = 1., int _sample_num = 160) {
        start = _start;
        goal = _goal;
        obstacle_list = obs;
        min_rand = rand_area[0];
        max_rand = rand_area[1];
        robot_radius = _robot_radius;
        sample_num = _sample_num;
        std::random_device seed;
        engine.seed(seed());
    }
    ~PRM() {}
    
    vector<vector<double>> planning(void);
};

vector<vector<double>> PRM::planning(void)
{
    vector<vector<double>> vertex = {{start[0], goal[0]}, {start[1], goal[1]}};

    if (show_animation) {
        plt::plot({start[0]}, {start[1]}, "xr");
        plt::plot({goal[0]}, {goal[1]}, "xr");
        for (vector<double> obs : obstacle_list) {
            plot_circle(obs[0], obs[1], obs[2], true);
        }
    }

    while (vertex[0].size() < sample_num + 2) {
        std::uniform_real_distribution<double> dist_xy(min_rand, max_rand);
        Vector2d rnd(dist_xy(engine), dist_xy(engine));
        if (check_collision(rnd)) {
            vertex[0].push_back(rnd[0]);
            vertex[1].push_back(rnd[1]);
        }
        if (show_animation && vertex[0].size() % 10 == 0) { 
            plt::plot(vertex[0], vertex[1], ".b");
            plt::pause(0.02);
        }
    }

    vector<vector<double>> edge(sample_num + 2, vector<double>(sample_num + 2, -1.));
    for (size_t i = 0; i < sample_num + 2; ++i) {
        for (size_t j = 0; j < sample_num + 2; ++j) {
            double node_distance = hypot(vertex[0][i] - vertex[0][j], vertex[1][i] - vertex[1][j]);
            if (node_distance < 3 &&
                check_collision({vertex[0][i], vertex[1][i]}, {vertex[0][j], vertex[1][j]})) {
                edge[i][j] = node_distance;
            }
        }
    }

    if (show_animation) {
        plot_road_map(vertex, edge);
    }

    vector<vector<double>> path = calc_final_path(vertex, edge);

    return path;
}

// Use dijkstra to implement single source shortest path search
vector<vector<double>> PRM::calc_final_path(
    const vector<vector<double>>& vertex, const vector<vector<double>>& edge)
{
    bool path_found = true;
    vector<vector<double>> path(2);
    Node start_node(vertex[0][0], vertex[1][0], 1000, -1);
    std::unordered_map<int, Node> open_set;
    std::unordered_map<int, Node> closed_set;
    open_set[0] = start_node;

    while (true) {
        if (open_set.size() == 0) {
            fmt::print("Cannot find path..\n");
            path_found = false;
            break;
        }

        int c_id = -1;
        for (const auto& n : open_set) {
            if (c_id == -1 || open_set[c_id].cost > n.second.cost) {
                c_id = n.first;
            }
        }
        Node current = open_set[c_id];
        open_set.erase(c_id);
        closed_set[c_id] = current;

        if (c_id == 1) {
            break;
        }
        
        for (int idx = 0; idx < edge.size(); ++idx) {
            if (edge[c_id][idx] > 0 && idx != c_id) {
                Node node(vertex[0][idx], vertex[1][idx],
                        current.cost + edge[c_id][idx], c_id);
                if (closed_set.find(idx) != closed_set.end()) {
                    continue;
                }
                if (open_set.find(idx) != open_set.end()) {
                    if (open_set[idx].cost > node.cost) {
                        open_set[idx].cost = node.cost;
                        open_set[idx].parent_index = c_id;
                    }
                } else {
                    open_set[idx] = node;
                }
            }
        }
    }

    if (!path_found) {
        return path;
    }

    int parent_index = 1;
    while (parent_index != -1) {
        Node n = closed_set[parent_index];
        path[0].push_back(n.x);
        path[1].push_back(n.y);
        parent_index = n.parent_index;
    }

    return path;
}

bool PRM::check_collision(Vector2d node)
{
    for (vector<double> obs : obstacle_list) {
        double dx = obs[0] - node[0];
        double dy = obs[1] - node[1];
        double d = hypot(dx, dy);

        if (d <= obs[2] + robot_radius) {
            return false;
        }
    }

    return true;
}

bool PRM::check_collision(Vector2d node1, Vector2d node2)
{
    if (hypot(node2[0] - node1[0], node2[1] - node1[1]) < 0.01) {
        return true;
    }

    for (vector<double> obs : obstacle_list) {
        double d1 = hypot(node1[0] - obs[0], node1[1] - obs[1]);
        double d2 = hypot(node2[0] - obs[0], node2[1] - obs[1]);

        if (d1 > d2) {
            std::swap(d1, d2);
            std::swap(node1, node2);
        }

        if (obs[2] >= d1 && obs[2] <= d2) {
            return false;
        } else if (obs[2] <= d1) {
            double d = abs( (node2[0] - node1[0]) * (node1[1] - obs[1]) -
                        (node1[0] - obs[0]) * (node2[1] - node1[1])) /
                        hypot(node2[0] - node1[0], node2[1] - node1[1]);
            Vector2d v1(obs[0] - node1[0], obs[1] - node1[1]);
            Vector2d v2(node2[0] - node1[0], node2[1] - node1[1]);
            if (d <= obs[2] && (v1[0] * v2[0] + v1[1] * v2[1]) >= 0) {
                return false;
            }
        }
    }

    return true;
}

void PRM::plot_road_map(
    const vector<vector<double>>& vertex, const vector<vector<double>>& edge)
{
    for (size_t i = 0; i < edge.size(); ++i) {
        for (size_t j = i + 1; j < edge.size(); ++j) {
            if (edge[i][j] > 0) {
                plt::plot({vertex[0][i], vertex[0][j]}, {vertex[1][i], vertex[1][j]}, "-g");
            }
        }
        if (i % 5 == 0) {
            plt::pause(0.02);
        }
    }
    // Avoid edge lines covering vertex
    plt::plot({start[0]}, {start[1]}, "xr");
    plt::plot({goal[0]}, {goal[1]}, "xr");
    plt::plot(vertex[0], vertex[1], ".b");
}

void PRM::plot_circle(double x, double y, double size, bool is_fill, string style)
{
    vector<double> xl;
    vector<double> yl;
    for (double deg = 0; deg <= M_PI * 2; deg += ( M_PI / 36.0)) {
        xl.push_back(x + size * cos(deg));
        yl.push_back(y + size * sin(deg));
    }
    if (is_fill) {
        plt::fill(xl, yl, {{"color", "gray"}});
    } else {
        plt::plot(xl, yl, style);
    }
}

int main(int argc, char** argv)
{
    vector<vector<double>> obstacle_list = {{5, 5, 1}, {3, 6, 2}, {3, 8, 2}, 
                            {3, 10, 2}, {7, 5, 2}, {9, 5, 2}, {8, 10, 1}};
    utils::TicToc t_m;
    Vector2d start(0, 0);
    Vector2d goal(6, 10);
    Vector2d area(-2, 13);
    vector<vector<double>> boundary(2);
    for (double i = area[0] - 1; i < (area[1] + 1); i += 0.2) {
        boundary[0].push_back(i);
        boundary[0].push_back(area[1] + 1);
        boundary[0].push_back(i);
        boundary[0].push_back(area[0] - 1);
        boundary[1].push_back(area[0] - 1);
        boundary[1].push_back(i);
        boundary[1].push_back(area[1] + 1);
        boundary[1].push_back(i);
    }

    if (show_animation) {
        plt::plot(boundary[0], boundary[1], "sk");
        plt::axis("equal");
        plt::grid(true);
        plt::title("Probabilistic RoadMap");
    }

    PRM prm(start, goal, obstacle_list, area);
    vector<vector<double>> path = prm.planning();
    if (show_animation) {
        plt::plot(path[0], path[1], "-r");
        plt::show();
    }
    
    return 0;
}
