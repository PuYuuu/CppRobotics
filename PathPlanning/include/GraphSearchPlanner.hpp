#pragma once
#ifndef __GRAPHSEARCHPLANNER_HPP
#define __GRAPHSEARCHPLANNER_HPP

#include <fmt/core.h>

#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

#include "utils.hpp"

class Node {
public:
    double x;
    double y;
    double cost;
    double parent_index;
    std::shared_ptr<Node> parent;

    Node() {}
    Node(double _x, double _y, double _cost = 0, int _parent_index = -1,
         std::shared_ptr<Node> _node = nullptr) {
        x = _x;
        y = _y;
        cost = _cost;
        parent_index = _parent_index;
        parent = _node;
    }
    ~Node() {}
};

class GraphSearchPlanner {
private:
    double minx;
    double miny;
    double maxx;
    double maxy;
    double xwidth;
    double ywidth;
    double map_resolution;
    double robot_radius;
    std::vector<std::vector<bool>> obstacle_map;
    std::vector<std::vector<double>> motion;

public:
    GraphSearchPlanner() {}
    GraphSearchPlanner(std::vector<double> ox, std::vector<double> oy, double reso, double radius) {
        map_resolution = reso;
        robot_radius = radius;
        calc_obstacle_map(ox, oy);
        motion = get_motion_model();
    }
    virtual ~GraphSearchPlanner() {}
    void calc_obstacle_map(const std::vector<double>& ox, const std::vector<double>& oy);
    std::vector<std::vector<double>> get_motion_model(void);
    double calc_grid_position(int index, double minp);
    double calc_xyindex(double position, double min_pos);
    double calc_grid_index(std::shared_ptr<Node> node);
    bool verify_node(std::shared_ptr<Node> node);
    std::vector<std::vector<double>> calc_final_path(
        std::shared_ptr<Node> ngoal, std::unordered_map<double, std::shared_ptr<Node>>& closed_set);

    virtual std::vector<std::vector<double>> planning(double sx, double sy, double gx,
                                                      double gy) = 0;

    double get_minx(void) const { return minx; }

    double get_miny(void) const { return miny; }

    std::vector<std::vector<double>> get_motion(void) const { return motion; }
};

#endif
