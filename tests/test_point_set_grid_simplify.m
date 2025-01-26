% test point_set_grid_simplify

clc;

addpath(genpath('../src'));
addpath('../data');


load('spiky_cell_like_surface.mat');
% load('unit_ball_16384_pts.mat'); % grid_step = 0.2

h = figure;
subplot(121);
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

marker = '+';
linewidth = 2;

scatter3(P(:,1),P(:,2),P(:,3),ones(size(P,1),1),C,marker,'LineWidth',linewidth);
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1],'FontSize',16);
view(2);
title('Raw point set','Color',[1 1 1],'FontSize',16);

grid_step = 0.1;
[G,C] = point_set_grid_simplify(P,grid_step,'exact',C);

subplot(122);
scatter3(G(:,1),G(:,2),G(:,3),ones(size(G,1),1),C,marker,'LineWidth',linewidth);
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1],'FontSize',16);
view(2);
title('Grid simplified point set','Color',[1 1 1],'FontSize',16);