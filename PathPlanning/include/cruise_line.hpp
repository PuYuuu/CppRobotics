#pragma once
#ifndef __CRUISE_LINE_HPP
#define __CRUISE_LINE_HPP

#include <cmath>
#include <vector>

class CruiseLine
{
private:
    double max_c;
    double road_width;

public:
    std::vector<std::vector<double>> ref_line;
    std::vector<std::vector<double>> bound_in;
    std::vector<std::vector<double>> bound_out;

    CruiseLine(double c = 0.15, double width = 8.0) :
        max_c(c), road_width(width) {}
    ~CruiseLine() {}

    std::vector<std::vector<double>> design_reference_line(void);
    std::vector<std::vector<double>> design_boundary_in(void);
    std::vector<std::vector<double>> design_boundary_out(void);
};

std::vector<std::vector<double>> CruiseLine::design_reference_line(void)
{
    std::vector<std::vector<double>> rxy(2);
    double step_curve = M_PI * 0.1;
    size_t step_line = 4;

    size_t cx = 30;
    size_t cy = 30;
    size_t cr = 20;
    for (double itheta = M_PI; itheta < 1.5 * M_PI; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double ix = 30; ix < 80; ix += step_line) {
        rxy[0].push_back(ix);
        rxy[1].push_back(10);
    }

    cx = 80;
    cy = 25;
    cr = 15;
    for (double itheta = -M_PI_2; itheta < M_PI_2; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double ix = 80; ix > 60; ix -= step_line) {
        rxy[0].push_back(ix);
        rxy[1].push_back(40);
    }

    cx = 60;
    cy = 60;
    cr = 20;
    for (double itheta = -M_PI_2; itheta > -M_PI; itheta -= step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }

    cx = 25;
    cy = 60;
    cr = 15;
    for (double itheta = 0.0; itheta < M_PI; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double iy = 60; iy > 30; iy -= step_line) {
        rxy[0].push_back(10);
        rxy[1].push_back(iy);
    }

    // rxy[0].push_back(rxy[0][0]);
    // rxy[1].push_back(rxy[1][0]);

    return rxy;
}

std::vector<std::vector<double>> CruiseLine::design_boundary_in(void)
{
    std::vector<std::vector<double>> rxy(2);
    double step_curve = M_PI * 0.1;
    size_t step_line = 2;

    size_t cx = 30;
    size_t cy = 30;
    size_t cr = 20 - road_width;
    for (double itheta = M_PI; itheta < 1.5 * M_PI; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double ix = 30; ix < 80; ix += step_line) {
        rxy[0].push_back(ix);
        rxy[1].push_back(10 + road_width);
    }

    cx = 80;
    cy = 25;
    cr = 15 - road_width;
    for (double itheta = -M_PI_2; itheta < M_PI_2; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double ix = 80; ix > 60; ix -= step_line) {
        rxy[0].push_back(ix);
        rxy[1].push_back(40 - road_width);
    }

    cx = 60;
    cy = 60;
    cr = 20 + road_width;
    for (double itheta = -M_PI_2; itheta > -M_PI; itheta -= step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }

    cx = 25;
    cy = 60;
    cr = 15 - road_width;
    for (double itheta = 0.0; itheta < M_PI; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double iy = 60; iy > 30; iy -= step_line) {
        rxy[0].push_back(10 + road_width);
        rxy[1].push_back(iy);
    }

    rxy[0].push_back(rxy[0][0]);
    rxy[1].push_back(rxy[1][0]);

    return rxy;
}

std::vector<std::vector<double>> CruiseLine::design_boundary_out(void)
{
    std::vector<std::vector<double>> rxy(2);
    double step_curve = M_PI * 0.05;
    size_t step_line = 2;

    size_t cx = 30;
    size_t cy = 30;
    size_t cr = 20 + road_width;
    for (double itheta = M_PI; itheta < 1.5 * M_PI; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double ix = 30; ix < 80; ix += step_line) {
        rxy[0].push_back(ix);
        rxy[1].push_back(10 - road_width);
    }

    cx = 80;
    cy = 25;
    cr = 15 + road_width;
    for (double itheta = -M_PI_2; itheta < M_PI_2; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double ix = 80; ix > 60; ix -= step_line) {
        rxy[0].push_back(ix);
        rxy[1].push_back(40 + road_width);
    }

    cx = 60;
    cy = 60;
    cr = 20 - road_width;
    for (double itheta = -M_PI_2; itheta > -M_PI; itheta -= step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }

    cx = 25;
    cy = 60;
    cr = 15 + road_width;
    for (double itheta = 0.0; itheta < M_PI; itheta += step_curve) {
        rxy[0].emplace_back(cx + cr * cos(itheta));
        rxy[1].emplace_back(cy + cr * sin(itheta));
    }
    for (double iy = 60; iy > 30; iy -= step_line) {
        rxy[0].push_back(10 - road_width);
        rxy[1].push_back(iy);
    }

    rxy[0].push_back(rxy[0][0]);
    rxy[1].push_back(rxy[1][0]);

    return rxy;
}

#endif
