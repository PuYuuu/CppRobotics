#include "dynamic_programming_heuristic.hpp"

#include <cmath>
#include <queue>
#include <unordered_map>

#include "utils.hpp"

using std::unordered_map;
using std::vector;

class NNode {
public:
    int x;
    int y;
    double cost;
    int pind;

    NNode(int _x, int _y, double _c, int _p) : x(_x), y(_y), cost(_c), pind(_p) {}
    NNode() {}
    ~NNode() {}
    bool operator<(const NNode& other) const { return cost > other.cost; }
};

class NPara {
public:
    int minx;
    int miny;
    int maxx;
    int maxy;
    int xw;
    int yw;
    double reso;
    vector<vector<double>> motion;

    NPara(int _minx, int _miny, int _maxx, int _maxy, int _xw, int _yw, double _r,
          vector<vector<double>> _m)
        : minx(_minx),
          miny(_miny),
          maxx(_maxx),
          maxy(_maxy),
          xw(_xw),
          yw(_yw),
          reso(_r),
          motion(_m) {}
    ~NPara() {}
};

vector<vector<double>> get_motion(void) {
    vector<vector<double>> motion = {{-1, 0}, {-1, 1}, {0, 1},  {1, 1},
                                     {1, 0},  {1, -1}, {0, -1}, {-1, -1}};

    return motion;
}

int calc_index(NNode node, const NPara& P) { return (node.y - P.miny) * P.xw + (node.x - P.minx); }

double u_cost(vector<double> u) { return hypot(u[0], u[1]); }

vector<vector<int>> calc_obsmap(const vector<double>& ox, const vector<double>& oy, double rr,
                                const NPara& P) {
    vector<vector<int>> obsmap(P.xw, vector<int>(P.yw, 0));
    for (size_t x = 0; x < P.xw; ++x) {
        int xx = x + P.minx;
        for (size_t y = 0; y < P.yw; ++y) {
            int yy = y + P.miny;
            for (size_t idx = 0; idx < ox.size(); ++idx) {
                if (hypot(ox[idx] - xx, oy[idx] - yy) <= rr / P.reso) {
                    obsmap[x][y] = 1;
                    break;
                }
            }
        }
    }

    return obsmap;
}

NPara calc_parameters(const vector<double>& ox, const vector<double>& oy, double rr, double reso,
                      vector<vector<int>>& obsmap) {
    int minx = round(utils::min(ox));
    int miny = round(utils::min(oy));
    int maxx = round(utils::max(ox));
    int maxy = round(utils::max(oy));
    int xw = maxx - minx;
    int yw = maxy - miny;

    vector<vector<double>> motion = get_motion();
    NPara P(minx, miny, maxx, maxy, xw, yw, reso, motion);
    obsmap = calc_obsmap(ox, oy, rr, P);

    return P;
}

bool check_node(NNode n, NPara P, const vector<vector<int>>& obsmap) {
    if (n.x <= P.minx || n.x >= P.maxx || n.y <= P.miny || n.y >= P.maxy ||
        obsmap[n.x - P.minx][n.y - P.miny]) {
        return false;
    }

    return true;
}

vector<vector<double>> calc_holonomic_heuristic_with_obstacle(std::shared_ptr<Node> node,
                                                              const vector<vector<double>>& obs,
                                                              double reso, double rr) {
    NNode n_goal(round(node->x.back() / reso), round(node->y.back() / reso), 0.0, -1);
    vector<double> ox;
    vector<double> oy;
    for (size_t idx = 0; idx < obs[0].size(); ++idx) {
        ox.push_back(obs[0][idx] / reso);
        oy.push_back(obs[1][idx] / reso);
    }
    vector<vector<int>> obsmap;
    NPara P = calc_parameters(ox, oy, reso, rr, obsmap);
    unordered_map<int, NNode> open_set;
    unordered_map<int, NNode> closed_set;
    open_set[calc_index(n_goal, P)] = n_goal;

    std::priority_queue<NNode> q_priority;
    q_priority.push(n_goal);

    while (true) {
        if (open_set.empty()) {
            break;
        }
        NNode n_curr = q_priority.top();
        q_priority.pop();
        int ind = calc_index(n_curr, P);
        closed_set[ind] = n_curr;
        open_set.erase(ind);

        for (int i = 0; i < P.motion.size(); ++i) {
            NNode nnode(n_curr.x + P.motion[i][0], n_curr.y + P.motion[i][1],
                        n_curr.cost + u_cost(P.motion[i]), ind);

            if (!check_node(nnode, P, obsmap)) {
                continue;
            }

            int n_ind = calc_index(nnode, P);
            if (closed_set.find(n_ind) == closed_set.end()) {
                if (open_set.find(n_ind) != open_set.end()) {
                    if (open_set[n_ind].cost > nnode.cost) {
                        open_set[n_ind].cost = nnode.cost;
                        open_set[n_ind].pind = ind;
                    }
                } else {
                    open_set[n_ind] = nnode;
                    q_priority.push(nnode);
                }
            }
        }
    }

    vector<vector<double>> hmap(P.xw, vector<double>(P.yw, std::numeric_limits<double>::max()));

    for (const std::pair<int, NNode>& kv : closed_set) {
        NNode n = kv.second;
        hmap[n.x - P.minx][n.y - P.miny] = n.cost;
    }

    return hmap;
}
