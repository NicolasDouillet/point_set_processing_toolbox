% test sinusoidal_icosahedron

clc;

addpath(genpath('../src'));
addpath('../data');

nb_samples = 120;
P = sinusoidal_icosahedron(nb_samples);
plot_point_set(P);