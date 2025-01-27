% test sinusoidal_dodecahedron

clc;

addpath(genpath('../src'));
addpath('../data');

nb_samples = 120;
w = 3; % shape parameter

P = sinusoidal_dodecahedron(nb_samples,w);
plot_point_set(P);