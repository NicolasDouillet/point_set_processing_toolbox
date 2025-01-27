% test sinusoidal_icosahedron

clc;

addpath(genpath('../src'));
addpath('../data');

V = sinusoidal_icosahedron(120);
plot_point_set(V);