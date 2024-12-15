% test grid_simplify

clc;

addpath(genpath('../src'));
addpath('../data');

load('virus_cell_like_surface.mat');
% load('unit_ball_16384_pts.mat'); % grid_step = 0.2


h = figure;
subplot(121);
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

marker = '.';
color = 'y';

plot3(P(:,1),P(:,2),P(:,3),marker,'Color',color,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(2);
title('Raw point set','Color',[1 1 1],'FontSize',16);


grid_step = 0.1;
S = point_set_grid_simplify(P,grid_step);


subplot(122);
plot3(S(:,1),S(:,2),S(:,3),marker,'Color',color,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(2);
title('Grid simplified point set','Color',[1 1 1],'FontSize',16);