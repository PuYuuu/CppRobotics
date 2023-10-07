#pragma once

#ifndef __GRAPHSEARCHPLANNER_HPP
#define __GRAPHSEARCHPLANNER_HPP

#include <cmath>
#include <vector>
#include <memory>
#include <unordered_map>

#include <fmt/core.h>

#include "utils/utils.hpp"

using std::vector;
using std::shared_ptr;
using std::unordered_map;

class Node
{
public:
    double x;
    double y;
    double cost;
    double parent_index;
    shared_ptr<Node> parent; 

    Node() {}
    Node(double _x, double _y, double _cost = 0, 
        int _parent_index = -1, shared_ptr<Node> _node = nullptr) {
        x = _x;
        y = _y;
        cost = _cost;
        parent_index = _parent_index;
        parent = _node;
    }
    ~Node() {}
};

class GraphSearchPlanner
{
private:
    double minx;
    double miny;
    double maxx;
    double maxy;
    double xwidth;
    double ywidth;
    double map_resolution;
    double robot_radius;
    vector<vector<bool>> obstacle_map;
    vector<vector<double>> motion;
public:
    GraphSearchPlanner() {}
    GraphSearchPlanner(vector<double> ox, vector<double> oy, double reso, double radius) {
        map_resolution = reso;
        robot_radius = radius;
        calc_obstacle_map(ox, oy);
        motion = get_motion_model();
    }
    ~GraphSearchPlanner() {}
    void calc_obstacle_map(const vector<double>& ox, const vector<double>& oy);
    vector<vector<double>> get_motion_model(void);
    double calc_grid_position(int index, double minp);
    double calc_xyindex(double position, double min_pos);
    double calc_grid_index(shared_ptr<Node> node);
    bool verify_node(shared_ptr<Node> node);
    vector<vector<double>> calc_final_path(shared_ptr<Node> ngoal, 
        unordered_map<double, shared_ptr<Node>>& closed_set);
    
    virtual vector<vector<double>> planning(double sx, double sy, double gx, double gy) = 0;

    double get_minx(void) const {
        return minx;
    }
    double get_miny(void) const {
        return miny;
    }
    vector<vector<double>> get_motion(void) const {
        return motion;
    }
};

#endif
