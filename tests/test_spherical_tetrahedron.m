% test spherical_tetrahedron

clc;

addpath(genpath('../src'));
addpath('../data');


Rho = 8;
r = 1.5;
nb_samples = 180;

P = spherical_tetrahedron(Rho,r,nb_samples);
plot_point_set(P);