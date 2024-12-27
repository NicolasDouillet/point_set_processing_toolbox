% test point_set_random_simplify

clc;

addpath(genpath('../src'));
addpath('../data');


load('unit_ball_4096_pts.mat');

plot_point_set(P,C);
title('Raw point set');
coeff = 0.5;
[P,C] = point_set_random_simplify(P,coeff,C);
plot_point_set(P,C);
title('Random simplified point set');