% test sinusoidal_dodecahedron

clc;

addpath(genpath('../src'));
addpath('../data');

V = sinusoidal_dodecahedron(120,3);
plot_point_set(V);