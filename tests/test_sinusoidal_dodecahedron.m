% test sinusoidal_dodecahedron

clc;

addpath(genpath('../src'));
addpath('../data');

nb_samples = 64;
w = 3; % shape parameter
random_sampling = false;

P = sinusoidal_dodecahedron(nb_samples,w,random_sampling);
plot_point_set(P);