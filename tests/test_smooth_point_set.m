% test smooth_point_set

clc;

addpath(genpath('../src'));
addpath('../data');


load('noisy_kitten_5%.mat');
% load('polar_10%_noisy_sphere_4096_pts.mat');
% load('unit_sphere_4096_pts_noisy.mat');
marker = '.';
color = 'y';
markersize = 6;
    
h = figure;
subplot(131);
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

plot3(P(:,1),P(:,2),P(:,3),marker,'Color',color,'MarkerSize',markersize,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);
title('Noisy kitten','Color',[1 1 1],'FontSize',16);

k = 51;
nb_it = 2;
P = smooth_point_set(P,k,nb_it);

subplot(132);
plot3(P(:,1),P(:,2),P(:,3),marker,'Color',color,'MarkerSize',markersize,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);
title('Smoothed kitten','Color',[1 1 1],'FontSize',16);


load('kitten.mat');


subplot(133);
plot3(P(:,1),P(:,2),P(:,3),marker,'Color',color,'MarkerSize',markersize,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);
title('Original kitten','Color',[1 1 1],'FontSize',16);