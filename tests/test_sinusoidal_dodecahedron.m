% test sinusoidal_dodecahedron

clc;

addpath(genpath('../src'));
addpath('../data');

nb_samples = 180;
w = 3; % shape parameter
random_sampling = true;

P = sinusoidal_dodecahedron(nb_samples,w,random_sampling);
plot_point_set(P);