#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <queue>
#include <unordered_map>

#include <fmt/core.h>

#include "utils/utils.hpp"
#include "utils/matplotlibcpp.h"
#include "GraphSearchPlanner.hpp"

using std::pair;
using std::make_pair;
using std::vector;
using std::shared_ptr;
using std::unordered_map;
using namespace Eigen;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

struct cmp{
   bool operator()(pair<Node, Vector2d>& a, pair<Node, Vector2d>& b){
        if (a.second[0] == b.second[0]) {
            return a.second[1] > b.second[1];
        }

       return a.second[0] > b.second[0]; 
   }
};

class DStarLite
{
private:
    int x_min_world;
    int y_min_world;
    int x_max;
    int y_max;
    Node start;
    Node goal;
    double km;
    Vector2d kold;
    bool initialized;
    vector<vector<double>> rhs;
    vector<vector<double>> g;
    vector<Node> obstacles;
    vector<Node> detected_obstacles;
    vector<double> detected_obstacles_for_plotting_x;
    vector<double> detected_obstacles_for_plotting_y;
    vector<vector<double>> motion;
    std::priority_queue<pair<Node, Vector2d>, vector<pair<Node, Vector2d>>, cmp> U;
public:
    DStarLite(const vector<double>& ox, const vector<double>& oy);
    ~DStarLite() {}
    vector<vector<double>> create_grid(double val);
    void run(Node st, Node go, const vector<double>& spoofed_ox, 
            const vector<double>& spoofed_oy);
    void initialize(Node st, Node go);
    Vector2d calculate_key(Node s);
    double h(Node s) {
        // return std::max(abs(start.x - s.x), abs(start.y - s.y));
        return 1.0;
    }
    void compute_shortest_path(void);
    bool compare_keys(Vector2d key_pair1, Vector2d key_pair2);
    void update_vertex(Node u);
    bool compare_coordinates(Node node1, Node node2);
    double c(Node node1, Node node2);
};

DStarLite::DStarLite(const vector<double>& ox, const vector<double>& oy)
{
    x_min_world = static_cast<int>(Utils::min(ox));
    y_min_world = static_cast<int>(Utils::min(oy));
    x_max = static_cast<int>(abs(Utils::max(ox) - x_min_world));
    y_max = static_cast<int>(abs(Utils::max(oy) - y_min_world));
    for (int idx = 0; idx < ox.size(); ++idx) {
        Node obs(ox[idx], oy[idx]);
        obstacles.emplace_back(obs);
    }
    start = Node(0, 0);
    goal = Node(0, 0);
    km = 0.0;
    kold << 0.0, 0.0;
    rhs = create_grid(std::numeric_limits<double>::max() / 2);
    g = create_grid(std::numeric_limits<double>::max() / 2);
    vector<vector<double>> motion_tmp = {{1, 0, 1}, {0, 1, 1}, {-1, 0, 1}, {0, -1, 1},
                {-1, -1, sqrt(2)}, {-1, 1, sqrt(2)}, {1, -1, sqrt(2)}, {1, 1, sqrt(2)}};
    motion = motion_tmp;
    initialized = false;
}

vector<vector<double>> DStarLite::create_grid(double val)
{
    vector<vector<double>> grid(x_max, vector<double>(y_max, val));
    return grid;
}

void DStarLite::run(Node st, Node go, const vector<double>& spoofed_ox,
                    const vector<double>& spoofed_oy)
{
    vector<double> pathx;
    vector<double> pathy;
    initialize(st, go);
    Node last = start;
    compute_shortest_path();

    pathx.push_back(start.x + x_min_world);
    pathy.push_back(start.y + y_min_world);
    if (show_animation) {
        // current_path = compute_current_path();
        // previous_path = current_path.copy()
        // previous_path_image = display_path(previous_path, ".c", alpha=0.3)
        // current_path_image = display_path(current_path, ".c")
    }
}   

void DStarLite::initialize(Node st, Node go)
{
    start.x = st.x - x_min_world;
    start.y = st.y - y_min_world;
    goal.x = go.x - x_min_world;
    goal.y = go.y - y_min_world;

    if (!initialized) {
        initialized = true;
        fmt::print("Initializing...");
        rhs[goal.x][goal.y] = 0;
        U.push(make_pair(goal, calculate_key(goal)));
    }
}

Vector2d DStarLite::calculate_key(Node s)
{
    Vector2d key;
    key << std::min(g[s.x][s.y], rhs[s.x][s.y]) + h(s) + km,
           std::min(g[s.x][s.y], rhs[s.x][s.y]);
    
    return key;
}

void DStarLite::compute_shortest_path(void)
{
    bool has_elements = (U.size() > 0);
    bool start_key_not_updated = compare_keys(
            U.top().second, calculate_key(start));
    bool rhs_not_equal_to_g = rhs[start.x][start.y] != g[start.x][start.y];
    
    while (has_elements && (start_key_not_updated || rhs_not_equal_to_g)) {
        kold = U.top().second;
        Node u = U.top().first;
        U.pop();
        if (compare_keys(kold, calculate_key(u))) {
            U.push(make_pair(u, calculate_key(u)));
        } else if (g[u.x][u.y] > rhs[u.x][u.y]) {
            g[u.x][u.y] = rhs[u.x][u.y];
            for (const vector<double>& m : motion) {
                Node node(u.x + m[0], u.y + m[1], u.cost + m[2]);
                if (0 <= node.x && node.x < x_max && 0 <= node.y && node.y < y_max) {
                    update_vertex(node);
                }
            }
        } else {
            g[u.x][u.y] = std::numeric_limits<double>::max();
            update_vertex(u);
            for (const vector<double>& m : motion) {
                Node node(u.x + m[0], u.y + m[1], u.cost + m[2]);
                if (0 <= node.x && node.x < x_max && 0 <= node.y && node.y < y_max) {
                    update_vertex(node);
                }
            }
        }
        start_key_not_updated = compare_keys(U.top().second, calculate_key(start));
        rhs_not_equal_to_g = rhs[start.x][start.y] != g[start.x][start.y];
    }
}

bool DStarLite::compare_keys(Vector2d key_pair1, Vector2d key_pair2)
{
    return (key_pair1[0] < key_pair2[0] ||
            (key_pair1[0] == key_pair2[0] && key_pair1[1] < key_pair2[1]));
}

void DStarLite::update_vertex(Node u)
{
    if (!compare_coordinates(u, goal)) {
        for (const vector<double>& m : motion) {
            Node node(u.x + m[0], u.y + m[1], u.cost + m[2]);
            if (0 <= node.x && node.x < x_max && 0 <= node.y && node.y < y_max) {
                rhs[u.x][u.y] = std::min(rhs[u.x][u.y], g[u.x][u.y] + c(u, node));
            }

            vector<pair<Node, Vector2d>> container;
            while (!U.empty()) {
                pair<Node, Vector2d> tmp = U.top();
                U.pop();
                if (!compare_coordinates(u, tmp.first)) {
                    container.emplace_back(tmp);
                } else {
                    break;
                }
            }
            for (auto con : container) {
                U.push(con);
            }

            if (g[u.x][u.y] != rhs[u.x][u.y]) {
                U.push(make_pair(u, calculate_key(u)));
            }
        }
    }
}

bool DStarLite::compare_coordinates(Node node1, Node node2)
{
    return node1.x == node2.x && node1.y == node2.y;
}

double DStarLite::c(Node node1, Node node2)
{

}

int main(int argc, char** argv)
{
    Node start(10, 10);
    Node goal(50, 50);
    vector<double> obstacle_x;
    vector<double> obstacle_y;
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
        plt::plot({start.x}, {start.y}, "og");
        plt::plot({goal.x}, {goal.y}, "xb");
        plt::grid(true);
        plt::axis("equal");

        plt::plot();
        plt::pause(0.01);
    }

    vector<double> spoofed_ox;
    vector<double> spoofed_oy;
    for (int i = 0; i < 21; ++i) {
        spoofed_ox.emplace_back(i);
        spoofed_oy.emplace_back(20);
    }
    for (int i = 0; i < 20; ++i) {
        spoofed_ox.emplace_back(0);
        spoofed_oy.emplace_back(i);
    }
    
    DStarLite dstarlite(obstacle_x, obstacle_y);

    return 0;
}
