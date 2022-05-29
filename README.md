# Bspline Trajectory Server with Optimization (bs-traj-server)

## Installation
`bs-traj-server` serves a complete trajectory server that uses bspline to represent the path, several test nodes are compiled using CMake.
- Includes a RRT-on-point-clouds as a front-end search, for an initial path
- Unconstrained Optimization (LBFGS-B) to acquire new control points

### Dependencies
- libbspline (https://github.com/matthewoots/libbspline.git) 
- LBFGSpp (https://github.com/matthewoots/LBFGSpp.git)

### Setup
```bash
git clone https://github.com/matthewoots/bs-traj-server.git --recurse-submodules
cd <bs-traj-server directory>
mkdir build && cd build
cmake .. 
make
```

#### Run executable
Run `./bs_traj_server_test_random_point_node` in the `build` folder, this is to test the performance and speed of and simulating test points and commands given
- *2 methods are tested here, a method to use idx and another is by getting from the time difference (faster approach)*
- *This can only be launched after you create and compile your **build** folder*
- *Output example is as shown below*
```bash
[tserver] test chronos: 0.00107656s
[tserver] cp_size/acceptable 63/25
...
[tserver] path_size: 180
test (1) runtime: 0.015661s
displaced_time: 2.0149
[tserver] cp_size/acceptable 63/25
...
[tserver] path_size: 180
[tserver] cp_size/acceptable 63/25
...
update_get_command_on_path_by_idx: 
 -1.11535
-0.460259
   1.7103
cmd_time: 2.04444
cmd_difference: 0.0114077
update_get_command_by_time: 
 -1.12347
-0.460554
  1.72431
cmd_time: 2.03254
cmd_difference: 4.611e-06
```

