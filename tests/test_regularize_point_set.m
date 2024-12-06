% test regularize_point_set

clc;

addpath(genpath('../src'));
addpath('../data');

    
load('unit_sphere_4096_vtx_noisy.mat');


disp('Radius variance before regularization : ');
var(sqrt(sum(V.^2,2)))
plot_point_set(V), view(2);

nb_ngb = 15;
nb_it = 3;
V = regularize_point_set(V,nb_ngb,nb_it); 

disp('Radius variance after regularization : ');
var(sqrt(sum(V.^2,2)))
plot_point_set(V), view(2);