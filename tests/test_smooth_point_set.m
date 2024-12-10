% test smooth_point_set

clc;

addpath(genpath('../src'));
addpath('../data');


load('unit_sphere_4096_vtx_noisy.mat');
marker = '.';
color = 'y';
markersize = 6;
    
h = figure;
subplot(121);
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

plot3(V(:,1),V(:,2),V(:,3),marker,'Color',color,'MarkerSize',markersize,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(2);
title('Noisy sphere','Color',[1 1 1],'FontSize',16);

k = 51;
nb_it = 2;
V = smooth_point_set(V,k,nb_it);

subplot(122);
plot3(V(:,1),V(:,2),V(:,3),marker,'Color',color,'MarkerSize',markersize,'MarkerEdgeColor',color,'MarkerFaceColor',color), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(2);
title('Smoothed sphere','Color',[1 1 1],'FontSize',16);