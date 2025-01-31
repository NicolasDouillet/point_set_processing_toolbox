% test convex Reuleaux tetrahedron

clc;

addpath(genpath('../src'));
addpath('../data');

nb_samples = 32;
random_sampling = false;

P = Reuleaux_tetrahedron(nb_samples,random_sampling);
plot_point_set(P);