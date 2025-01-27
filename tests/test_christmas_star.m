% test christmas star

clc;

addpath(genpath('../src'));
addpath('../data');

V = christmas_star(32);
plot_point_set(V);