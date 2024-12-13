% test regularize_point_set

clc;

addpath(genpath('../src'));
addpath('../data');

    
load('unit_sphere_4096_pts_noisy.mat');


disp('Radius variance before regularization : ');
var(sqrt(sum(P.^2,2)))
plot_point_set(P), view(2);

nb_ngb = 15;
nb_it = 3;
P = regularize_point_set(P,nb_ngb,nb_it); 

disp('Radius variance after regularization : ');
var(sqrt(sum(P.^2,2)))
plot_point_set(P), view(2);