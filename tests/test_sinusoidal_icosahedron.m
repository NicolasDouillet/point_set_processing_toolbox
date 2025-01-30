% test sinusoidal_icosahedron

clc;

addpath(genpath('../src'));
addpath('../data');

nb_samples = 32;
w = 1; % shape parameter
random_sampling = false;

P = sinusoidal_icosahedron(nb_samples,w,random_sampling);
plot_point_set(P);