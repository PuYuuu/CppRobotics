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
constexpr double N_STEER = 10.0;                // steer command number
constexpr double GOAL_YAW_ERROR = M_PI / 60;
constexpr double XY_RESO = 2.0;                 // [m]
constexpr double YAW_RESO = 15 * M_PI / 180;    // [rad]
constexpr double MOVE_STEP = 0.2;               // [m] path interporate resolution
constexpr double COLLISION_CHECK_STEP = 10;     // skip number for collision check
constexpr double EXTEND_AREA = 5.0;             // collision check range extended
constexpr double GEAR_COST = 100.0;             // switch back penalty cost
constexpr double BACKWARD_COST = 5.0;           // backward penalty cost
constexpr double STEER_CHANGE_COST = 5.0;       // steer angle change penalty cost
constexpr double STEER_ANGLE_COST = 1.0;        // steer angle penalty cost
constexpr double SCISSORS_COST = 200.0;         // scissors cost
constexpr double H_COST = 10.0;                 // Heuristic cost penalty cost

class Para
{
public:
    int minx;
    int miny;
    int minyaw;
    int minyawt;
    int maxx;
    int maxy;
    int maxyaw;
    int maxyawt;
    int xw;
    int yw;
    int yaww;
    int yawtw;
    double xyreso;
    double yawreso;
    vector<vector<double>> obs;
    KDTree* kdtree;
    utils::VehicleConfig vc;

    Para(int _minx, int _miny, int _minyaw, int _minyawt, int _maxx, int _maxy,
        int _maxyaw, int _maxyawt, int _xw, int _yw, int _yaww, int _yawtw, 
        double _xyreso, double _yawreso, vector<vector<double>> _obs, KDTree* _kdtree,
        utils::VehicleConfig _vc) :
        minx(_minx), miny(_miny), minyaw(_minyaw), minyawt(_minyawt), maxx(_maxx), 
        maxy(_maxy), maxyaw(_maxyaw), maxyawt(_maxyawt), xw(_xw), yw(_yw), yaww(_yaww),
        yawtw(_yawtw), xyreso(_xyreso), yawreso(_yawreso), obs(_obs), kdtree(_kdtree), vc(_vc) {

    }
    ~Para () {}
};

struct cmp{
    bool operator() ( pair<int, double> a, pair<int, double> b ) {
        return a.second > b.second;
    }
};

vector<vector<double>> generate_obstacle(void)
{
    vector<vector<double>> obs(2);
    for (int i = -30; i < 31; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(38);
    }
    for (int i = -30; i < -6; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(23);
    }
    for (int i = 7; i < 31; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(23);
    }
    for (int i = 0; i < 24; ++i) {
        obs[0].push_back(-6);
        obs[1].push_back(i);
    }
    for (int i = 0; i < 24; ++i) {
        obs[0].push_back(6);
        obs[1].push_back(i);
    }
    for (int i = -6; i < 7; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(0);
    }

    return obs;
}

vector<vector<double>> generate_obstacle_2(void)
{
    vector<vector<double>> obs(2);
    for (int i = -30; i < 31; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(20);
    }
    for (int i = -30; i < -11; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(8);
    }
    for (int i = 12; i < 31; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(8);
    }
    for (int i = 0; i < 8; ++i) {
        obs[0].push_back(-12);
        obs[1].push_back(i);
    }
    for (int i = 0; i < 8; ++i) {
        obs[0].push_back(12);
        obs[1].push_back(i);
    }
    for (int i = -12; i < 12; ++i) {
        obs[0].push_back(i);
        obs[1].push_back(0);
    }

    return obs;
}

Para calc_parameters(vector<vector<double>> obs,
    double xyreso, double yawreso, KDTree* kdtree, utils::VehicleConfig vc)
{
    double minxm = utils::min(obs[0]) - EXTEND_AREA;
    double minym = utils::min(obs[1]) - EXTEND_AREA;
    double maxxm = utils::max(obs[0]) + EXTEND_AREA;
    double maxym = utils::max(obs[1]) + EXTEND_AREA;

    obs[0].push_back(minxm);
    obs[1].push_back(minym);
    obs[0].push_back(maxxm);
    obs[1].push_back(maxym);

    int minx = round(minxm / xyreso);
    int miny = round(minym / xyreso);
    int maxx = round(maxxm / xyreso);
    int maxy = round(maxym / xyreso);

    int xw = maxx - minx;
    int yw = maxy - miny;
    int minyaw = round(-M_PI / yawreso) - 1;
    int maxyaw = round(M_PI / yawreso);
    int yaww = maxyaw - minyaw;

    return Para(minx, miny, minyaw, minyaw, maxx, maxy, maxyaw, maxyaw,
                xw, yw, yaww, yaww, xyreso, yawreso, obs, kdtree, vc);
}

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

int calc_index(const shared_ptr<const Node>& node, const Para& P)
{
    int ind = (node->yawind - P.minyaw) * P.xw * P.yw +
              (node->yind - P.miny) * P.xw + (node->xind - P.minx);

    int yawt_ind = round(node->yawt.back() / P.yawreso);
    ind += (yawt_ind - P.minyawt) * P.xw * P.yw * P.yaww;

    return ind;
}

double calc_hybrid_cost(
    const shared_ptr<const Node>& node, const vector<vector<double>>& hmap, const Para& P)
{
    double cost = node->cost + H_COST * hmap[node->xind - P.minx][node->yind - P.miny];

    return cost;
}

vector<double> calc_trailer_yaw(
    const vector<double>& yaw, double yawt0, const vector<double>& steps, const Para& P)
{
    vector<double> yawt(yaw.size(), 0.0);
    yawt[0] = yawt0;

    for (size_t idx = 1; idx < yaw.size(); ++idx) {
        yawt[idx] += 
            (yawt[idx - 1] + steps[idx - 1] / P.vc.RTR * sin(yaw[idx - 1] - yawt[idx - 1]));
    }

    return yawt;
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
    
    for (size_t idx = 0; idx < rspath.yaw.size(); ++idx) {
        cost += SCISSORS_COST * abs(utils::pi_2_pi(rspath.yaw[idx] - rspath.yawt[idx]));
    }

    return cost;
}

bool is_collision(vector<double>& x, vector<double>& y,
    vector<double>& yaw, vector<double>& yawt, const Para& P)
{
    for (size_t idx = 0; idx < x.size(); ++idx) {
        double d = 1;
        double deltal = (P.vc.RTF - P.vc.RTB) / 2.0;
        double rt = (P.vc.RTF + P.vc.RTB) / 2.0 + d;
        double ctx = x[idx] + deltal * cos(yaw[idx]);
        double cty = y[idx] + deltal * sin(yaw[idx]);
        vector<point_t> ids = P.kdtree->neighborhood_points({ctx, cty}, rt);

        if (ids.empty()) {
            continue;
        } 
        for (const point_t& ob : ids) {
            double xot = ob[0] - ctx;
            double yot = ob[1] - cty;
            double dx_trail = xot * cos(yawt[idx]) + yot * sin(yawt[idx]);
            double dy_trail = -xot * sin(yawt[idx]) + yot * cos(yawt[idx]);
            if (abs(dx_trail) <= rt && abs(dy_trail) <= P.vc.W / 2.0 + d) {
                return true;
            }
        }

        deltal = (P.vc.RF - P.vc.RB) / 2.0;
        double rc = (P.vc.RF + P.vc.RB) / 2.0 + d;
        double cx = x[idx] + deltal * cos(yaw[idx]);
        double cy = y[idx] + deltal * sin(yaw[idx]);
        ids = P.kdtree->neighborhood_points({cx, cy}, rc);

        if (ids.empty()) {
            continue;
        } 
        for (const point_t& ob : ids) {
            double xo = ob[0] - cx;
            double yo = ob[1] - cy;
            double dx_car = xo * cos(yaw[idx]) + yo * sin(yaw[idx]);
            double dy_car = -xo * sin(yaw[idx]) + yo * cos(yaw[idx]);

            if (abs(dx_car) <= rc && abs(dy_car) <= P.vc.W / 2 + d) {
                return true;
            }
        }
    }

    return false;
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
    for (Path& path : paths) {
        vector<double> steps;
        for (int d : path.directions) {
            steps.push_back(MOVE_STEP * d);
        }
        vector<double> yawt = calc_trailer_yaw(path.yaw, node->yawt.back(), steps, P);
        path.yawt = yawt;
    }
    std::sort(paths.begin(), paths.end(), [&P](const Path& a, const Path& b) {
        return calc_rs_path_cost(a, P) < calc_rs_path_cost(b, P);
    });

    Path path = paths[0];
    vector<double> pathx;
    vector<double> pathy;
    vector<double> pathyaw;
    vector<double> pathyawt;
    for (size_t idx = 0; idx < path.x.size(); idx += COLLISION_CHECK_STEP) {
        pathx.push_back(path.x[idx]);
        pathy.push_back(path.y[idx]);
        pathyaw.push_back(path.yaw[idx]);
        pathyawt.push_back(path.yawt[idx]);
    }
    if (!is_collision(pathx, pathy, pathyaw, pathyawt, P)) {
        return path;
    }

    return Path();
}

bool update_node_with_analystic_expantion(shared_ptr<Node> n_curr, 
    shared_ptr<Node> ngoal, double gyawt, shared_ptr<Node>& fpath, const Para& P)
{
    Path path = analystic_expantion(n_curr, ngoal, P);
    if (path.x.empty() || path.y.empty() || path.yaw.empty()) {
        return false;
    }

    if (abs(utils::pi_2_pi(path.yawt.back() - gyawt)) >= GOAL_YAW_ERROR) {
        return false;
    }

    vector<double> fx(path.x.begin() + 1, path.x.end() - 1);
    vector<double> fy(path.y.begin() + 1, path.y.end() - 1);
    vector<double> fyaw(path.yaw.begin() + 1, path.yaw.end() - 1);
    vector<double> fyawt(path.yawt.begin() + 1, path.yawt.end() - 1);
    vector<int> fd(path.directions.begin() + 1, path.directions.end() - 1);
    double fcost = n_curr->cost + calc_rs_path_cost(path, P);
    int fpind = calc_index(n_curr, P);
    double fsteer = 0.0;
    fpath = make_shared<Node>(n_curr->xind, n_curr->yind, n_curr->yawind,
        n_curr->direction, fx, fy, fyaw, fyawt, fd, fsteer, fcost, fpind);

    return true;
}

bool is_index_ok(int xind, int yind, const vector<double>& xlist, 
    const vector<double>& ylist, const vector<double>& yawlist, 
    const vector<double>& yawtlist, const Para& P)
{
    if (xind <= P.minx || xind >= P.maxx || yind <= P.miny || yind >= P.maxy) {
        return false;
    }

    vector<double> nodex;
    vector<double> nodey;
    vector<double> nodeyaw;
    vector<double> nodeyawt;
    for (size_t idx = 0; idx < xlist.size(); idx += COLLISION_CHECK_STEP) {
        nodex.push_back(xlist[idx]);
        nodey.push_back(ylist[idx]);
        nodeyaw.push_back(yawlist[idx]);
        nodeyawt.push_back(yawtlist[idx]);
    }

    if (is_collision(nodex, nodey, nodeyaw, nodeyawt, P)) {
        return false;
    }

    return true;
}

shared_ptr<Node> calc_next_node(
    const shared_ptr<const Node>& n_curr, int c_id,double u, int d, const Para& P)
{
    double step = XY_RESO * 2.5;
    int nlist = ceil(step / MOVE_STEP);
    vector<double> xlist = {n_curr->x.back() + d * MOVE_STEP * cos(n_curr->yaw.back())};
    vector<double> ylist = {n_curr->y.back() + d * MOVE_STEP * sin(n_curr->yaw.back())};
    vector<double> yawlist = {utils::pi_2_pi(n_curr->yaw.back() + d * MOVE_STEP / P.vc.WB * tan(u))};
    vector<double> yawtlist = {utils::pi_2_pi(n_curr->yawt.back() +
                           d * MOVE_STEP / P.vc.RTR * sin(n_curr->yaw.back() - n_curr->yawt.back()))};
    shared_ptr<Node> node;

    for (size_t idx = 0; idx < nlist - 1; ++idx) {
        xlist.push_back(xlist[idx] + d * MOVE_STEP * cos(yawlist[idx]));
        ylist.push_back(ylist[idx] + d * MOVE_STEP * sin(yawlist[idx]));
        yawlist.push_back(utils::pi_2_pi(yawlist[idx] + d * MOVE_STEP / P.vc.WB * tan(u)));
        yawtlist.push_back(utils::pi_2_pi(yawtlist[idx] +
                            d * MOVE_STEP / P.vc.RTR * sin(yawlist[idx] - yawtlist[idx])));
    }
    int xind = round(xlist.back() / P.xyreso);
    int yind = round(ylist.back() / P.xyreso);
    int yawind = round(yawlist.back() / P.yawreso);

    if (!is_index_ok(xind, yind, xlist, ylist, yawlist, yawtlist, P)) {
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
    for (size_t idx = 0; idx < yawlist.size(); ++idx) {
        cost += SCISSORS_COST * abs(utils::pi_2_pi(yawlist[idx] - yawtlist[idx]));
    }

    cost = n_curr->cost + cost;
    vector<int> directions(xlist.size(), direction);
    node = make_shared<Node>(xind, yind, yawind, direction, xlist,
            ylist, yawlist, yawtlist, directions, u, cost, c_id);

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
    vector<double> ryawt;
    vector<int> direc;
    double cost = 0.0;
    shared_ptr<Node> node = ngoal;

    while (true) {
        rx.insert(rx.end(), node->x.rbegin(), node->x.rend());
        ry.insert(ry.end(), node->y.rbegin(), node->y.rend());
        ryaw.insert(ryaw.end(), node->yaw.rbegin(), node->yaw.rend());
        ryawt.insert(ryawt.end(), node->yawt.rbegin(), node->yawt.rend());
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
    std::reverse(ryawt.begin(), ryawt.end());
    std::reverse(direc.begin(), direc.end());
    direc[0] = direc[1];

    Path path(rx, ry, ryaw, ryawt, direc);

    return path;
}


Path hybrid_astar_planning(Vector4d start, Vector4d goal,
    vector<vector<double>>& obs, utils::VehicleConfig VC, double xyreso, double yawreso)
{
    int sxr = round(start[0] / xyreso);
    int syr = round(start[1] / xyreso);
    int syawr = round(utils::pi_2_pi(start[2]) / yawreso);
    int gxr = round(goal[0] / xyreso);
    int gyr = round(goal[1] / xyreso);
    int gyawr = round(utils::pi_2_pi(goal[2]) / yawreso);

    shared_ptr<Node> nstart(new Node(
        sxr, syr, syawr, 1, {start[0]}, {start[1]}, {start[2]}, {start[3]}, {1}, 0.0, 0.0, -1));
    shared_ptr<Node> ngoal(new Node(
        gxr, gyr, gyawr, 1, {goal[0]}, {goal[1]}, {goal[2]}, {goal[3]}, {1}, 0.0, 0.0, -1));
    
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

        bool update = update_node_with_analystic_expantion(n_curr, ngoal, goal[3], fnode, P);
        if (update) {
            break;
        }

        double yawt0 = n_curr->yawt[0];
        for (size_t idx = 0; idx < steer_set.size(); ++idx) {
            shared_ptr<Node> node = calc_next_node(n_curr, ind, steer_set[idx], direc_set[idx], P);

            if (node == nullptr) {
                continue;
            }

            // This will cause the program to take a very long time
            // if (show_animation) {
            //     plt::plot(node->x, node->y);
            //     plt::pause(0.001);
            // }

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
    Vector4d start;
    Vector4d goal;
    vector<vector<double>> obs;

    if (argc < 2) {
        fmt::print("Trailer truck parallel parking scene\n");
        start << -20, 14.0, 0, 0;
        goal << 2.5, 4, 0, 0;
        obs = generate_obstacle_2();
        plt::ylim(-15, 35);
        plt::title("Hybrid A* with Trailer - Case 0");
    } else {
        fmt::print("Trailer truck reversing into parking spot scene\n");
        start << 18.0, 34.0, M_PI, M_PI;
        goal << 0.0, 12.0, M_PI_2, M_PI_2;
        obs = generate_obstacle();
        plt::title("Hybrid A* with Trailer - Case 1");
    }
    utils::VehicleConfig vc(4.5, 1.0, 3.0, 3.5, 0.5, 1.0, 0.6);

    plt::plot(obs[0], obs[1], "sk");
    utils::draw_trailer(start, 0.0, vc);
    utils::draw_trailer(goal, 0.0, vc, "0.4");
    plt::pause(1);

    utils::TicToc t_m;
    Path path = hybrid_astar_planning(start, goal, obs, vc, XY_RESO, YAW_RESO);
    fmt::print("hybrid_astar_with_trailer planning costtime: {:.3f} s\n", t_m.toc() / 1000);

    if (path.x.empty() || path.y.empty() || path.yaw.empty()) {
        fmt::print("Searching failed!\n");
        return 0;
    }

    double steer = 0.0;
    for (size_t idx = 0; idx < path.x.size(); ++idx) {
        plt::cla();
        plt::plot(obs[0], obs[1], "sk");
        plt::plot(path.x, path.y, "-r");

        if (idx < path.x.size() - 2) {
            double dy = (path.yaw[idx + 1] - path.yaw[idx]) / MOVE_STEP;
            steer = utils::pi_2_pi(atan(vc.WB * dy / path.directions[idx]));
        } else {
            steer = 0.0;
        }

        utils::draw_trailer({path.x[idx], path.y[idx], path.yaw[idx], path.yawt[idx]}, steer, vc);
        if (idx < path.x.size() - 1) {
            utils::draw_trailer(goal, 0.0, vc, "0.4");
        }
        
        if (argc < 2) {
            plt::title("Hybrid A* with Trailer - Case 0");
        } else {
            plt::title("Hybrid A* with Trailer - Case 1");
        }
        plt::axis("equal");
        plt::pause(0.01);
    }
    plt::show();

    return 0;
}
