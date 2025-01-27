% test christmas star

clc;

addpath(genpath('../src'));
addpath('../data');


nb_samples = 32;
P = christmas_star(nb_samples);
plot_point_set(P);