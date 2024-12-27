% test plot_point_set (with or without colors)

clc;

addpath(genpath('../src'));
addpath('../data');


load('left_hand.mat');
C = C/255; % From |[0; 255]| to |[0; 1]| color scale
% plot_point_set(P); % without colors
plot_point_set(P,C);