#pragma once
#ifndef __DYNAMIC_PROGRAMMING_HEURISTIC_HPP
#define __DYNAMIC_PROGRAMMING_HEURISTIC_HPP

#include <vector>
#include <memory>

class Node
{
public:
    int xind;
    int yind;
    int yawind;
    int direction;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> yaw;
    std::vector<double> yawt;
    std::vector<int> directions;
    double steer;
    double cost;
    int pind;

    Node (int _xi, int _yi, int _yawi, int _d, std::vector<double> _x,
        std::vector<double> _y, std::vector<double> _yaw, std::vector<int> _ds,
        double _s, double _c, int _p) : xind(_xi), yind(_yi),yawind(_yawi),
        direction(_d), x(_x), y(_y), yaw(_yaw), directions(_ds),
        steer(_s), cost(_c), pind(_p) {}
    Node (int _xi, int _yi, int _yawi, int _d, std::vector<double> _x,
        std::vector<double> _y, std::vector<double> _yaw, std::vector<double> _yawt,
        std::vector<int> _ds, double _s, double _c, int _p) : 
        xind(_xi), yind(_yi),yawind(_yawi), direction(_d), x(_x), y(_y), yaw(_yaw),
        yawt(_yawt),directions(_ds), steer(_s), cost(_c), pind(_p) {}
    Node() {}
    ~Node() {}
};

std::vector<std::vector<double>> calc_holonomic_heuristic_with_obstacle(
    std::shared_ptr<Node> node, const std::vector<std::vector<double>>& obs, 
    double reso, double rr);

#endif
