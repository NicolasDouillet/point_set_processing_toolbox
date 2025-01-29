% test spherical_tetrahedron

clc;

addpath(genpath('../src'));
addpath('../data');


Rho = 8;
r = 1.5;
nb_samples = 180;
random_sampling = true;

P = spherical_tetrahedron(Rho,r,nb_samples,random_sampling);
plot_point_set(P);