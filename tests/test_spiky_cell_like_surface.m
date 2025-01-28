% test spiky_cell_like_surface

clc;

addpath(genpath('../src'));
addpath('../data');

Rho = 1;
alpha_l = pi/12;
nb_samples = 240;
n = 4;
option_random_sampling = false;

P = spiky_cell_like_surface(Rho,alpha_l,nb_samples,n,option_random_sampling);
plot_point_set(P);
view(2);