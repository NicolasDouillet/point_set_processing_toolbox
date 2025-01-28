% test n_holes_torus

clc;

addpath(genpath('../src'));
addpath('../data');

a = 1;
b = 1;
Rho = 8;
r = 2;
n = 2;
nb_samples = 120;
option_random_sampling = false;

P = n_holes_torus(a,b,Rho,r,n,nb_samples,option_random_sampling);
plot_point_set(P);