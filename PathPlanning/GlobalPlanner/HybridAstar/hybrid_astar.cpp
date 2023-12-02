#include <cmath>
#include <vector>
#include <queue>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include <fmt/core.h>

#include "utils/utils.hpp"
#include "utils/KDTree.hpp"
#include "utils/matplotlibcpp.h"
#include "PathPlanning/CurvesGenerator/reeds_shepp_path.hpp"
#include "dynamic_programming_heuristic.hpp"

using std::vector;
using std::pair;
using std::shared_ptr;
using std::unordered_map;
using std::make_pair;
using std::make_shared;
using namespace Eigen;
namespace plt = matplotlibcpp;

constexpr bool show_animation = true;
constexpr int N_STEER = 3;                      // steer command number
constexpr double XY_RESO = 2.0;                 // [m]
constexpr double YAW_RESO = 15 * M_PI / 180;    // [rad]
constexpr double MOVE_STEP = 0.4;               // [m] path interporate resolution
constexpr double COLLISION_CHECK_STEP = 5;      // skip number for collision check
constexpr double EXTEND_BOUND = 1;              // collision check range extended
constexpr double GEAR_COST = 100.0;             // switch back penalty cost
constexpr double BACKWARD_COST = 5.0;           // backward penalty cost
constexpr double STEER_CHANGE_COST = 5.0;       // steer angle change penalty cost
constexpr double STEER_ANGLE_COST = 1.0;        // steer angle penalty cost
constexpr double H_COST = 15.0;                 // Heuristic cost penalty cost

class Para
{
public:
    int minx;
    int miny;
    int minyaw;
    int maxx;
    int maxy;
    int maxyaw;
    int xw;
    int yw;
    int yaww;
    double xyreso;
    double yawreso;
    vector<vector<double>> obs;
    KDTree* kdtree;
    utils::VehicleConfig vc;

    Para(int _minx, int _miny, int _minyaw, int _maxx, int _maxy,
        int _maxyaw, int _xw, int _yw, int _yaww, double _xyreso,
        double _yawreso, vector<vector<double>> _obs, KDTree* _kdtree,
        utils::VehicleConfig _vc) :
        minx(_minx), miny(_miny), minyaw(_minyaw), maxx(_maxx), maxy(_maxy),
        maxyaw(_maxyaw), xw(_xw), yw(_yw), yaww(_yaww), xyreso(_xyreso),
        yawreso(_yawreso), obs(_obs), kdtree(_kdtree), vc(_vc) {

    }
    ~Para () {}
};

struct cmp{
    bool operator() ( pair<int, double> a, pair<int, double> b ) {
        return a.second > b.second;
    }
};

pair<vector<double>, vector<int>> calc_motion_set(const utils::VehicleConfig& C)
{
    double step = C.MAX_STEER / N_STEER;
    vector<double> steer;
    vector<int> direc;

    for (double i = -C.MAX_STEER; i <= C.MAX_STEER; i += step) {
        steer.push_back(i);
    }
    for (size_t i = 0; i < steer.size(); ++i) {
        direc.push_back(1);
    }
    for (size_t i = 0; i < steer.size(); ++i) {
        direc.push_back(-1);
    }
    for (double i = -C.MAX_STEER; i <= C.MAX_STEER; i += step) {
        steer.push_back(i);
    }

    return make_pair(steer, direc);
}

vector<vector<double>> generate_obstacle(double x, double y)
{
    vector<vector<double>> obs(2);

    for (size_t idx = 0; idx < x; ++idx) {
        obs[0].push_back(idx);
        obs[1].push_back(0);
    }
    for (size_t idx = 0; idx < x; ++idx) {
        obs[0].push_back(idx);
        obs[1].push_back(y - 1);
    }
    for (double idx = 0; idx < y - 0.5; idx += 0.5) {
        obs[0].push_back(0);
        obs[1].push_back(idx);
    }
    for (double idx = 0; idx < y - 0.5; idx += 0.5) {
        obs[0].push_back(x - 1);
        obs[1].push_back(idx);
    }
    for (size_t idx = 10; idx < 21; ++idx) {
        obs[0].push_back(idx);
        obs[1].push_back(15);
    }
    for (double idx = 0; idx < 15; idx += 0.5) {
        obs[0].push_back(20);
        obs[1].push_back(idx);
    }
    for (double idx = 15; idx < 30; idx += 0.5) {
        obs[0].push_back(30);
        obs[1].push_back(idx);
    }
    for (double idx = 0; idx < 16; idx += 0.5) {
        obs[0].push_back(40);
        obs[1].push_back(idx);
    }
    
    return obs;
}

Para calc_parameters(vector<vector<double>> obs,
    double xyreso, double yawreso, KDTree* kdtree, utils::VehicleConfig vc)
{
    int minx = round(utils::min(obs[0]) / xyreso);
    int miny = round(utils::min(obs[1]) / xyreso);
    int maxx = round(utils::max(obs[0]) / xyreso);
    int maxy = round(utils::max(obs[1]) / xyreso);

    int xw = maxx - minx;
    int yw = maxy - miny;
    int minyaw = round(-M_PI / yawreso) - 1;
    int maxyaw = round(M_PI / yawreso);
    int yaww = maxyaw - minyaw;

    return Para(minx, miny, minyaw, maxx, maxy, maxyaw, xw,
                yw, yaww, xyreso, yawreso, obs, kdtree, vc);
}

int calc_index(const shared_ptr<const Node>& node, Para P)
{
    int ind = (node->yawind - P.minyaw) * P.xw * P.yw +
              (node->yind - P.miny) * P.xw + (node->xind - P.minx);

    return ind;
}

double calc_hybrid_cost(
    const shared_ptr<const Node>& node, const vector<vector<double>>& hmap, const Para& P)
{
    double cost = node->cost + H_COST * hmap[node->xind - P.minx][node->yind - P.miny];

    return cost;
}

bool is_collision(vector<double>& x, vector<double>& y, vector<double>& yaw, const Para& P)
{
    for (size_t idx = 0; idx < x.size(); ++idx) {
        int d = 1;
        double dl = (P.vc.RF - P.vc.RB) / 2.0;
        double r = (P.vc.RF + P.vc.RB) / 2.0 + d;
        double cx = x[idx] + dl * cos(yaw[idx]);
        double cy = y[idx] + dl * sin(yaw[idx]);
        vector<point_t> ids = P.kdtree->neighborhood_points({cx, cy}, r);

        for (const point_t& ob : ids) {
            double xo = ob[0] - cx;
            double yo = ob[1] - cy;
            double dx = xo * cos(yaw[idx]) + yo * sin(yaw[idx]);
            double dy = -xo * sin(yaw[idx]) + yo * cos(yaw[idx]);

            if (abs(dx) < r && abs(dy) < P.vc.W / 2 + d) {
                return true;
            }
        }
    }

    return false;
}

double calc_rs_path_cost(const Path& rspath, const Para& P)
{
    double cost = 0.0;

    for (double lr : rspath.lengths) {
        if (lr >= 0) {
            cost += 1;
        } else {
            cost += abs(lr) * BACKWARD_COST;
        }
    }
    for (size_t idx = 0; idx < rspath.lengths.size() - 1; idx++) {
        if (rspath.lengths[idx] * rspath.lengths[idx + 1] < 0.0) {
            cost += GEAR_COST;
        }
    }
    for (char ctype : rspath.ctypes) {
        if (ctype != 'S') {
            cost += STEER_ANGLE_COST * abs(P.vc.MAX_STEER);
        }
    }

    vector<double> ulist(rspath.ctypes.size(), 0.0);
    for (size_t idx = 0; idx < rspath.ctypes.size(); ++idx) {
        if (rspath.ctypes[idx] == 'R') {
            ulist[idx] = -P.vc.MAX_STEER;
        } else if (rspath.ctypes[idx] == 'L') {
            ulist[idx] = P.vc.MAX_STEER;
        }
    }
    for (size_t idx = 0; idx < rspath.ctypes.size() - 1; ++idx) {
        cost += (STEER_CHANGE_COST * abs(ulist[idx + 1] - ulist[idx]));
    }

    return cost;
}

Path analystic_expantion(
    const shared_ptr<const Node>& node, const shared_ptr<const Node>& ngoal, const Para& P)
{
    Vector3d start(node->x.back(), node->y.back(), node->yaw.back());
    Vector3d goal(ngoal->x.back(), ngoal->y.back(), ngoal->yaw.back());
    double maxc = tan(P.vc.MAX_STEER) / P.vc.WB;
    vector<Path> paths = calc_rs_paths(start, goal, maxc, MOVE_STEP);
    
    if (paths.empty()) {
        return Path();
    }
    std::sort(paths.begin(), paths.end(), [&P](const Path& a, const Path& b) {
        return calc_rs_path_cost(a, P) < calc_rs_path_cost(b, P);
    });

    for (Path path : paths) {
        vector<double> pathx;
        vector<double> pathy;
        vector<double> pathyaw;
        for (size_t idx = 0; idx < path.x.size(); idx += COLLISION_CHECK_STEP) {
            pathx.push_back(path.x[idx]);
            pathy.push_back(path.y[idx]);
            pathyaw.push_back(path.yaw[idx]);
        }
        if (!is_collision(pathx, pathy, pathyaw, P)) {
            return path;
        }
    }

    return Path();
}

bool update_node_with_analystic_expantion(
    shared_ptr<Node> n_curr, shared_ptr<Node> ngoal, shared_ptr<Node>& fpath, const Para& P)
{
    Path path = analystic_expantion(n_curr, ngoal, P);

    if (path.x.empty() || path.y.empty() || path.yaw.empty()) {
        return false;
    }

    vector<double> fx(path.x.begin() + 1, path.x.end() - 1);
    vector<double> fy(path.y.begin() + 1, path.y.end() - 1);
    vector<double> fyaw(path.yaw.begin() + 1, path.yaw.end() - 1);
    vector<int> fd(path.directions.begin() + 1, path.directions.end() - 1);
    double fcost = n_curr->cost + calc_rs_path_cost(path, P);
    int fpind = calc_index(n_curr, P);
    double fsteer = 0.0;
    fpath = make_shared<Node>(n_curr->xind, n_curr->yind, n_curr->yawind,
        n_curr->direction, fx, fy, fyaw, fd, fsteer, fcost, fpind);

    return true;
}

bool is_index_ok(int xind, int yind, const vector<double>& xlist, 
    const vector<double>& ylist, const vector<double>& yawlist, const Para& P)
{
    if (xind <= P.minx || xind >= P.maxx || yind <= P.miny || yind >= P.maxy) {
        return false;
    }

    vector<double> nodex;
    vector<double> nodey;
    vector<double> nodeyaw;
    for (size_t idx = 0; idx < xlist.size(); idx += COLLISION_CHECK_STEP) {
        nodex.push_back(xlist[idx]);
        nodey.push_back(ylist[idx]);
        nodeyaw.push_back(yawlist[idx]);
    }

    if (is_collision(nodex, nodey, nodeyaw, P)) {
        return false;
    }

    return true;
}

shared_ptr<Node> calc_next_node(
    const shared_ptr<const Node>& n_curr, int c_id, double u, int d, const Para& P)
{
    double step = XY_RESO * 2.5;
    int nlist = ceil(step / MOVE_STEP);
    vector<double> xlist = {n_curr->x.back() + d * MOVE_STEP * cos(n_curr->yaw.back())};
    vector<double> ylist = {n_curr->y.back() + d * MOVE_STEP * sin(n_curr->yaw.back())};
    vector<double> yawlist = {utils::pi_2_pi(n_curr->yaw.back() + d * MOVE_STEP / P.vc.WB * tan(u))};
    shared_ptr<Node> node;

    for (size_t idx = 0; idx < nlist - 1; ++idx) {
        xlist.push_back(xlist[idx] + d * MOVE_STEP * cos(yawlist[idx]));
        ylist.push_back(ylist[idx] + d * MOVE_STEP * sin(yawlist[idx]));
        yawlist.push_back(utils::pi_2_pi(yawlist[idx] + d * MOVE_STEP / P.vc.WB * tan(u)));
    }
    int xind = round(xlist.back() / P.xyreso);
    int yind = round(ylist.back() / P.xyreso);
    int yawind = round(yawlist.back() / P.yawreso);

    if (!is_index_ok(xind, yind, xlist, ylist, yawlist, P)) {
        return nullptr;
    }

    double cost = 0.0;
    int direction = 1;
    if (d > 0) {
        direction = 1;
        cost += abs(step);
    } else {
        direction = -1;
        cost += abs(step) * BACKWARD_COST;
    }
    if (direction != n_curr->direction) {
        cost += GEAR_COST;
    }
    cost += STEER_ANGLE_COST * abs(u);
    cost += STEER_CHANGE_COST * abs(n_curr->steer - u);
    cost = n_curr->cost + cost;
    vector<int> directions(xlist.size(), direction);
    node = make_shared<Node>(xind, yind, yawind, direction, xlist,
            ylist, yawlist, directions, u, cost, c_id);

    return node;
}

bool is_same_grid(const shared_ptr<const Node>& node1, const shared_ptr<const Node>& node2)
{
    if (node1->xind != node2->xind || node1->yind != node2->yind || node1->yawind != node2->yawind) {
        return false;
    }

    return true;
}

Path extract_path(unordered_map<int, shared_ptr<Node>>& closed,
                shared_ptr<Node> ngoal, shared_ptr<Node> nstart)
{
    vector<double> rx;
    vector<double> ry;
    vector<double> ryaw;
    vector<int> direc;
    double cost = 0.0;
    shared_ptr<Node> node = ngoal;

    while (true) {
        rx.insert(rx.end(), node->x.rbegin(), node->x.rend());
        ry.insert(ry.end(), node->y.rbegin(), node->y.rend());
        ryaw.insert(ryaw.end(), node->yaw.rbegin(), node->yaw.rend());
        direc.insert(direc.end(), node->directions.rbegin(), node->directions.rend());
        cost += node->cost;

        if (is_same_grid(node, nstart)) {
            break;
        }

        node = closed[node->pind];
    }

    std::reverse(rx.begin(), rx.end());
    std::reverse(ry.begin(), ry.end());
    std::reverse(ryaw.begin(), ryaw.end());
    std::reverse(direc.begin(), direc.end());
    direc[0] = direc[1];

    Path path(rx, ry, ryaw, direc);

    return path;
}

Path hybrid_astar_planning(Vector3d start, Vector3d goal,
    const vector<vector<double>>& obs, utils::VehicleConfig VC, double xyreso, double yawreso)
{
    int sxr = round(start[0] / xyreso);
    int syr = round(start[1] / xyreso);
    int syawr = round(utils::pi_2_pi(start[2]) / yawreso);
    int gxr = round(goal[0] / xyreso);
    int gyr = round(goal[1] / xyreso);
    int gyawr = round(utils::pi_2_pi(goal[2]) / yawreso);

    shared_ptr<Node> nstart(new Node(
        sxr, syr, syawr, 1, {start[0]}, {start[1]}, {start[2]}, {1}, 0.0, 0.0, -1));
    shared_ptr<Node> ngoal(new Node(
        gxr, gyr, gyawr, 1, {goal[0]}, {goal[1]}, {goal[2]}, {1}, 0.0, 0.0, -1));
    pointVec points;
    for (size_t idx = 0; idx < obs[0].size(); ++idx) {
        points.push_back({obs[0][idx], obs[1][idx]});
    }
    KDTree kdtree(points);
    
    Para P = calc_parameters(obs, xyreso, yawreso, &kdtree, VC);

    vector<vector<double>> hmap =
        calc_holonomic_heuristic_with_obstacle(ngoal, P.obs, P.xyreso, 1.0);
    pair<vector<double>, vector<int>> steer_and_dir = calc_motion_set(P.vc);
    vector<double> steer_set = steer_and_dir.first;
    vector<int> direc_set = steer_and_dir.second;
    unordered_map<int, shared_ptr<Node>> open_set;
    unordered_map<int, shared_ptr<Node>> closed_set;
    open_set[calc_index(nstart, P)] = nstart;
    std::priority_queue<pair<int, double>, vector<pair<int, double>>, cmp> q_priority;
    q_priority.push(make_pair(calc_index(nstart, P), calc_hybrid_cost(nstart, hmap, P)));

    shared_ptr<Node> fnode = nullptr;
    while (true) {
        if (open_set.empty()) {
            return Path();
        }
        int ind = q_priority.top().first;
        q_priority.pop();
        shared_ptr<Node> n_curr = open_set[ind];
        closed_set[ind] = n_curr;
        open_set.erase(ind);

        bool update = update_node_with_analystic_expantion(n_curr, ngoal, fnode, P);
        if (update) {
            break;
        }

        for (size_t idx = 0; idx < steer_set.size(); ++idx) {
            shared_ptr<Node> node = calc_next_node(n_curr, ind, steer_set[idx], direc_set[idx], P);
            
            if (node == nullptr) {
                continue;
            }

            if (show_animation) {
                plt::plot(node->x, node->y, "-");
                plt::pause(0.001);
            }

            int node_ind = calc_index(node, P);
            if (closed_set.find(node_ind) != closed_set.end()) {
                continue;
            }

            if (open_set.find(node_ind) == open_set.end()) {
                open_set[node_ind] = node;
                q_priority.push(make_pair(node_ind, calc_hybrid_cost(node, hmap, P)));
            } else {
                if (open_set[node_ind]->cost > node->cost) {
                    open_set[node_ind] = node;
                }
            }
        }
    }
    fmt::print("final expand node: {}\n", open_set.size() + closed_set.size());
    
    return extract_path(closed_set, fnode, nstart);
}

int main(int argc, char** argv)
{
    Vector3d start(10.0, 7.0, 120 * M_PI / 180);
    Vector3d goal(45.0, 20.0, M_PI_2);
    vector<vector<double>> obs = generate_obstacle(51, 31);
    utils::VehicleConfig VC(4.5, 1.0, 3.0, 3.5, 0.5, 1.0, 0.6);

    plt::plot(obs[0], obs[1], "sk");
    utils::draw_vehicle(start, 0.0, VC);
    utils::draw_vehicle(goal, 0.0, VC, "0.4");
    plt::title("Hybrid A*");
    plt::pause(1.0);

    utils::TicToc t_m;
    Path path = hybrid_astar_planning(start, goal, obs, VC, XY_RESO, YAW_RESO);
    fmt::print("hybrid_astar planning costtime: {:.3f} s\n", t_m.toc() / 1000);

    if (path.x.empty() || path.y.empty() || path.yaw.empty()) {
        fmt::print("Searching failed!\n");
        return 0;
    }

    double steer = 0.0;
    for (size_t idx = 0; idx < path.x.size(); ++idx) {
        plt::cla();
        plt::plot(obs[0], obs[1], "sk");
        plt::plot(path.x, path.y, "r");

        if (idx < path.x.size() - 2) {
            double dy = (path.yaw[idx + 1] - path.yaw[idx]) / MOVE_STEP;
            steer = -utils::pi_2_pi(atan(-3.5 * dy / path.directions[idx]));
        } else {
            steer = 0.0;
        }

        utils::draw_vehicle({path.x[idx], path.y[idx], path.yaw[idx]}, steer, VC);
        if (idx < path.x.size() - 1) {
            utils::draw_vehicle(goal, 0.0, VC, "0.4");
        }
        plt::title("Hybrid A*");
        plt::axis("equal");
        plt::pause(0.01);
    }

    plt::cla();
    plt::plot(obs[0], obs[1], "sk");
    utils::draw_vehicle(start, 0.0, VC, "0.4");
    utils::draw_vehicle(goal, 0.0, VC);
    plt::title("Hybrid A*");
    plt::axis("equal");
    for (size_t idx = 2; idx < path.x.size(); idx += 4) {
        utils::draw_vehicle({path.x[idx], path.y[idx], path.yaw[idx]}, 0, VC, "-c", false, false);
        plt::pause(0.001);
    }

    plt::show();

    return 0;
}
