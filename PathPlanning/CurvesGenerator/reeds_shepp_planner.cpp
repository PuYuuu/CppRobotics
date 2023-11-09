#include <cmath>
#include <string>
#include <fmt/core.h>

#include "utils/utils.hpp"
#include "reeds_shepp_path.hpp"
#include "utils/matplotlibcpp.h"

using std::string;
namespace plt = matplotlibcpp;
constexpr bool show_animation = true;

int main(int argc, char** argv)
{
    Vector3d start(-10.0, -10.0, M_PI_4);
    Vector3d goal(0., 0., -M_PI_2);
    double curvature = 0.1;
    double step_size = 0.05;
    utils::VehicleConfig vc(0.5);

    Path path = reeds_shepp_path(start, goal, curvature, step_size);
    string final_mode = "final course ";
    final_mode.push_back(path.ctypes[0]);
    final_mode.push_back(path.ctypes[1]);
    final_mode.push_back(path.ctypes[2]);

    if (show_animation) {
        for (size_t idx = 0; idx < path.x.size(); ++idx) {
            plt::cla();
            plt::named_plot(final_mode, path.x, path.y);
            plt::arrow(start[0], start[1], cos(start[2]), sin(start[2]), "r", 0.075);
            plt::arrow(goal[0], goal[1], cos(goal[2]), sin(goal[2]), "g", 0.075);

            utils::draw_vehicle({path.x[idx], path.y[idx], path.yaw[idx]}, 0, vc, false);
            plt::legend({{"loc", "upper left"}});
            plt::grid(true);
            plt::axis("equal");
            plt::title("Reeds Shepp Path Planning");
            plt::pause(0.001);
        }
        plt::show();
    }

    return 0;
}
