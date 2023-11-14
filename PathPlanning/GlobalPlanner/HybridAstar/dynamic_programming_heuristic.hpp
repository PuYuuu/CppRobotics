#pragma once
#ifndef DYNAMIC_PROGRAMMING_HEURISTIC_HPP
#define DYNAMIC_PROGRAMMING_HEURISTIC_HPP

#include <vector>
#include <memory>

using std::vector;

class Node
{
public:
    int xind;
    int yind;
    int yawind;
    int direction;
    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    vector<int> directions;
    double steer;
    double cost;
    int pind;

    Node (int _xi, int _yi, int _yawi, int _d, vector<double> _x,
        vector<double> _y, vector<double> _yaw, vector<int> _ds,
        double _s, double _c, int _p) : xind(_xi), yind(_yi),yawind(_yawi),
        direction(_d), x(_x), y(_y), yaw(_yaw), directions(_ds),
        steer(_s), cost(_c), pind(_p) {

    }
    Node() {}
    ~Node() {}
};

vector<vector<double>> calc_holonomic_heuristic_with_obstacle(
    std::shared_ptr<Node> node, const vector<vector<double>>& obs, 
    double reso, double rr);

#endif
