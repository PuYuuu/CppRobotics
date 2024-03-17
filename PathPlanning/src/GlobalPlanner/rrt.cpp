#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "matplotlibcpp.h"
#include "utils.hpp"

using std::string;
using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

class Node {
public:
    double x;
    double y;
    Node* parent;

    Node() : x(0), y(0), parent(nullptr) {}
    Node(double _x, double _y, Node* p = nullptr) : x(_x), y(_y), parent(p) {}
    ~Node() {}
};

class RRT {
private:
    Node* start;
    Node* goal;
    double min_rand;
    double max_rand;
    double expand_dis;
    double goal_sample_rate;
    int max_iter;
    double robot_radius;
    vector<Node*> node_list;
    vector<vector<double>> obstacle_list;
    vector<vector<double>> boundary;
    std::mt19937 engine;

    Node get_random_node(void);
    Node* get_nearest_node(Node rnd_node);
    Node* steer(Node* from_node, Node* to_node);
    Vector2d calc_distance_and_angle(Node* from_node, Node* to_node);
    bool check_collision(Node* new_node);
    void plot_circle(double x, double y, double size, bool is_fill = false, string style = "-b");
    void draw_graph(Node rnd);
    vector<vector<double>> generate_final_course(void);

public:
    RRT(Vector2d _start, Vector2d _goal, vector<vector<double>> obs, Vector2d rand_area,
        double expand = 0.5, double goal_sample = 0.5, int _max_iter = 1000,
        double _robot_radius = 0.5) {
        start = new Node(_start[0], _start[1]);
        goal = new Node(_goal[0], _goal[1]);
        min_rand = rand_area[0];
        max_rand = rand_area[1];
        expand_dis = expand;
        goal_sample_rate = goal_sample;
        max_iter = _max_iter;
        robot_radius = _robot_radius;
        obstacle_list = obs;
        std::random_device seed;
        engine.seed(seed());
    }
    ~RRT();

    vector<vector<double>> planning(void);
};

RRT::~RRT() {
    for (Node* n : node_list) {
        delete n;
    }
}

vector<vector<double>> RRT::planning(void) {
    boundary.resize(2);
    for (double i = min_rand - 1; i < (max_rand + 1); i += 0.2) {
        boundary[0].push_back(i);
        boundary[0].push_back(max_rand + 1);
        boundary[0].push_back(i);
        boundary[0].push_back(min_rand - 1);
        boundary[1].push_back(min_rand - 1);
        boundary[1].push_back(i);
        boundary[1].push_back(max_rand + 1);
        boundary[1].push_back(i);
    }

    node_list.push_back(start);

    for (int iter = 0; iter < max_iter; ++iter) {
        Node rnd_node = get_random_node();
        Node* nearest_node = get_nearest_node(rnd_node);
        Node* new_node = steer(nearest_node, &rnd_node);

        if (check_collision(new_node)) {
            node_list.emplace_back(new_node);
        } else {
            delete new_node;
        }

        if (show_animation && (iter % 5 == 0)) {
            draw_graph(rnd_node);
        }

        Vector2d d_angle = calc_distance_and_angle(node_list.back(), goal);
        if (d_angle[0] <= expand_dis) {
            goal->parent = node_list.back();
            break;
        }
    }

    return generate_final_course();
}

Node RRT::get_random_node(void) {
    Node rnd;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    if (dist(engine) > goal_sample_rate) {
        std::uniform_real_distribution<double> dist_xy(min_rand, max_rand);
        rnd = Node(dist_xy(engine), dist_xy(engine));
    } else {
        rnd = Node(goal->x, goal->y);
    }

    return rnd;
}

Node* RRT::get_nearest_node(Node rnd_node) {
    Node* min_node = new Node(node_list[0]->x, node_list[0]->y);
    double distance = hypot(min_node->x - rnd_node.x, min_node->y - rnd_node.y);
    for (int idx = 1; idx < node_list.size(); ++idx) {
        double tmp = hypot(node_list[idx]->x - rnd_node.x, node_list[idx]->y - rnd_node.y);
        if (tmp < distance) {
            min_node = node_list[idx];
            distance = tmp;
        }
    }

    return min_node;
}

Node* RRT::steer(Node* from_node, Node* to_node) {
    Node new_node(from_node->x, from_node->x);
    Vector2d d_angle = calc_distance_and_angle(from_node, to_node);
    double dist = std::min(expand_dis, d_angle[0]);
    Node* node_new =
        new Node(from_node->x + dist * cos(d_angle[1]), from_node->y + dist * sin(d_angle[1]));
    node_new->parent = from_node;

    return node_new;
}

Vector2d RRT::calc_distance_and_angle(Node* from_node, Node* to_node) {
    Vector2d d_angle;
    double dx = to_node->x - from_node->x;
    double dy = to_node->y - from_node->y;
    double d = hypot(dx, dy);
    double angle = atan2(dy, dx);
    d_angle << d, angle;

    return d_angle;
}

bool RRT::check_collision(Node* new_node) {
    if (new_node == nullptr) {
        return true;
    }

    for (vector<double> obs : obstacle_list) {
        double dx = obs[0] - new_node->x;
        double dy = obs[1] - new_node->y;
        double d = hypot(dx, dy);

        if (d <= (obs[2] + robot_radius)) {
            return false;
        }
    }

    return true;
}

void RRT::plot_circle(double x, double y, double size, bool is_fill, string style) {
    vector<double> xl;
    vector<double> yl;
    for (double deg = 0; deg <= M_PI * 2; deg += (M_PI / 36.0)) {
        xl.push_back(x + size * cos(deg));
        yl.push_back(y + size * sin(deg));
    }
    if (is_fill) {
        plt::fill(xl, yl, {{"color", "gray"}});
    } else {
        plt::plot(xl, yl, style);
    }
}

void RRT::draw_graph(Node rnd) {
    plt::clf();

    plt::plot(boundary[0], boundary[1], "sk");
    plt::plot({rnd.x}, {rnd.y}, "^k");
    if (robot_radius > 0.0) {
        plot_circle(rnd.x, rnd.y, robot_radius, false, "-r");
    }

    for (Node* n : node_list) {
        if (n->parent != nullptr) {
            plt::plot({n->x, n->parent->x}, {n->y, n->parent->y}, "-g");
        }
    }

    for (vector<double> obs : obstacle_list) {
        plot_circle(obs[0], obs[1], obs[2], true);
    }

    plt::plot({start->x}, {start->y}, "xr");
    plt::plot({goal->x}, {goal->y}, "xr");
    plt::axis("equal");
    plt::grid(true);
    plt::title("Rapid-exploration Random Tree");
    plt::pause(0.01);
}

vector<vector<double>> RRT::generate_final_course(void) {
    vector<vector<double>> final_path(2);
    Node* n = goal;

    while (n != nullptr) {
        final_path[0].emplace_back(n->x);
        final_path[1].emplace_back(n->y);
        n = n->parent;
    }

    return final_path;
}

int main(int argc, char** argv) {
    vector<vector<double>> obstacle_list = {{5, 5, 1}, {3, 6, 2}, {3, 8, 2}, {3, 10, 2},
                                            {7, 5, 2}, {9, 5, 2}, {8, 10, 1}};
    utils::TicToc t_m;
    Vector2d start(0, 0);
    Vector2d goal(6, 10);
    Vector2d area(-2, 13);
    RRT rrt(start, goal, obstacle_list, area);
    vector<vector<double>> path = rrt.planning();
    fmt::print("rrt planning costtime: {:.3f} s\n", t_m.toc() / 1000);

    if (path[0].size() < 2) {
        fmt::print("planning failed!\n");
        return 0;
    }

    if (show_animation) {
        plt::plot(path[0], path[1], "-r");
        plt::grid(true);
        plt::show();
    }

    return 0;
}
