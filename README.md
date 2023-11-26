# CppRobotics
![Language](https://img.shields.io/badge/language-c++-brightgreen)  ![Static Badge](https://img.shields.io/badge/eigen-3.3.7-blue)

Cpp implementation of robotics algorithms including localization, path planning, path tracking and control, inspired by [PythonRobotics](https://github.com/AtsushiSakai/PythonRobotics) and [MotionPlanning](https://github.com/zhm-real/MotionPlanning).

## 📌Requirement

- cmake
- Eigen3
- [fmt](https://github.com/fmtlib/fmt)

## 🛠Build

```shell
git clone git@github.com:PuYuuu/CppRobotics.git
cd CppRobotics
mkdir build && cd build
cmake ..
make -j6
```

Find all the executable files in **$workspace/bin**.

## 🎈Animations

### Perception

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/ekf.gif" width="400"/></a></td>
    <td><img src="./assets/images/rectangle_fitting.gif" width="400"/></a></td>
  </tr>
</table>
</div>

### PathPlanning

#### GlobalPlanner

##### Search_based Planning

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/bfs.gif" width="400"/></a></td>
    <td><img src="./assets/images/dfs.gif" width="400"/></a></td>
  </tr>
</table>
</div>
<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/dijkstra.gif" width="400"/></a></td>
    <td><img src="./assets/images/astar.gif" width="400"/></a></td>
  </tr>
</table>
</div>
<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/astar-bidirectional.gif" width="400"/></a></td>
    <td><img src="./assets/images/hybrid_astar.gif" width="400"/></a></td>
  </tr>
</table>
</div>
<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/hybrid_astar_trailer_0.gif" width="400"/></a></td>
    <td><img src="./assets/images/hybrid_astar_trailer.gif" width="400"/></a></td>
  </tr>
</table>
</div>

##### Sampling_based Planning

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/rrt.gif" width="400"/></a></td>
    <td><img src="./assets/images/probabilistic_road_map.gif" width="400"/></a></td>
  </tr>
</table>
</div>

#### LocalPlanner

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/dwa.gif" width="400"/></a></td>
    <td><img src="./assets/images/potential_field.gif" width="400"/></a></td>
  </tr>
</table>
</div>

#### CurvesGenerator

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/cubic.png" width="400"/></a></td>
    <td><img src="./assets/images/bspline.png" width="400"/></a></td>
  </tr>
</table>
</div>
<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/dubins.gif" width="400"/></a></td>
    <td><img src="./assets/images/bezier.gif" width="400"/></a></td>
  </tr>
</table>
</div>
<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/reeds_shepp.gif" width="400"/></a></td>
  </tr>
</table>
</div>

### PathTracking

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/pure_pursuit.gif" width="400"/></a></td>
    <td><img src="./assets/images/stanley.gif" width="400"/></a></td>
  </tr>
</table>
</div>

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/lqr-xy.gif" width="400"/></a></td>
  </tr>
</table>
</div>

### Control

<div align=center>
<table>
  <tr>
    <td><img src="./assets/images/move_to_pose.gif" width="400"/></a></td>
    <td><img src="./assets/images/move_to_pose_robots.gif" width="400"/></a></td>
  </tr>
</table>
</div>

## 🧾Licence

MIT
