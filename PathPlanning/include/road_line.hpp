#pragma once
#ifndef __ROAD_LINE_HPP
#define __ROAD_LINE_HPP

#include <cmath>
#include <vector>

class RoadLine {
protected:
    double road_width;

public:
    std::vector<std::vector<double>> ref_line;
    std::vector<std::vector<double>> bound_in;
    std::vector<std::vector<double>> bound_out;

    explicit RoadLine(double width = 8.0) : road_width(width) {}
    ~RoadLine() {}

    virtual std::vector<std::vector<double>> design_reference_line(void) = 0;
    virtual std::vector<std::vector<double>> design_boundary_left(void) = 0;
    virtual std::vector<std::vector<double>> design_boundary_right(void) = 0;
};

class CruiseRoadLine : public RoadLine {
private:
    double max_c;

public:
    explicit CruiseRoadLine(double c = 0.15, double width = 8.0) : max_c(c), RoadLine(width) {}
    ~CruiseRoadLine() {}

    std::vector<std::vector<double>> design_reference_line(void) override;
    std::vector<std::vector<double>> design_boundary_left(void) override;
    std::vector<std::vector<double>> design_boundary_right(void) override;
};

std::vector<std::vector<double>> CruiseRoadLine::design_reference_line(void) {
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

std::vector<std::vector<double>> CruiseRoadLine::design_boundary_left(void) {
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

std::vector<std::vector<double>> CruiseRoadLine::design_boundary_right(void) {
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

class StopRoadLine : public RoadLine {
public:
    explicit StopRoadLine(double width = 8.0) : RoadLine(width) {}
    ~StopRoadLine() {}

    std::vector<std::vector<double>> design_reference_line(void) override;
    std::vector<std::vector<double>> design_boundary_left(void) override;
    std::vector<std::vector<double>> design_boundary_right(void) override;
};

std::vector<std::vector<double>> StopRoadLine::design_reference_line(void) {
    std::vector<std::vector<double>> rxy(2);

    for (double i = 0.0; i < 60.0; i += 1.0) {
        rxy[0].push_back(i);
        rxy[1].push_back(0.0);
    }

    return rxy;
}

std::vector<std::vector<double>> StopRoadLine::design_boundary_left(void) {
    std::vector<std::vector<double>> rxy(2);

    for (double i = 0.0; i < 60.0; i += 1.0) {
        rxy[0].push_back(i);
        rxy[1].push_back(road_width);
    }

    return rxy;
}

std::vector<std::vector<double>> StopRoadLine::design_boundary_right(void) {
    std::vector<std::vector<double>> rxy(2);

    for (double i = 0.0; i < 60.0; i += 1.0) {
        rxy[0].push_back(i);
        rxy[1].push_back(-road_width);
    }

    return rxy;
}

#endif
