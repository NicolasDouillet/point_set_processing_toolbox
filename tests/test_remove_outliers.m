% test remove_outliers

clc;

addpath(genpath('../src'));
addpath('../data');

load('unit_sphere_4096_with_outliers.mat');
plot_point_set(V);

V = remove_outliers(V,20,1);
plot_point_set(V);