% test biscuit

clc;

addpath(genpath('../src'));
addpath('../data');

L = 16;
l = 9;
e = 5;
nb_samples = 32;
isotropic_sampling = true;
random_sampling = false;

P = biscuit(L,l,e,nb_samples,isotropic_sampling,random_sampling);
plot_point_set(P);