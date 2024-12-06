% test grid_simplify

clc;

addpath(genpath('../src'));
addpath('../data');
addpath('../point_set_grid_simplify');

% load('unit_ball_16384_vtx.mat');
load('virus_cell_like_surface.mat');

grid_step = 0.05;
S = point_set_grid_simplify(V,grid_step);

plot_point_set(V);
view(2);

size(S)

plot_point_set(S);
view(2);