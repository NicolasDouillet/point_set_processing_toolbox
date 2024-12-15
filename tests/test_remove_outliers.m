% test remove_outliers

clc;

addpath(genpath('../src'));
addpath('../data');

load('unit_sphere_4096_pts_with_outliers.mat');

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
view(3);
title('Unit sphere 4096 points with outlier groups','Color',[1 1 1],'FontSize',16);

P = remove_outliers(P,31,2.5);

subplot(122);
plot3(P(:,1),P(:,2),P(:,3),marker,'Color',color,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);
title('Unit sphere 4096 points onces outliers removed','Color',[1 1 1],'FontSize',16);