% test biscuit

clc;

addpath(genpath('../src'));
addpath('../data');

L = 16;
l = 9;
e = 5;
nb_samples = 64;
random_sampling = false;

P = biscuit(L,l,e,nb_samples,random_sampling);
plot_point_set(P);
view(2);