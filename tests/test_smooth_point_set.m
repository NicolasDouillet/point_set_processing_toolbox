% test smooth_point_set

clc;

addpath(genpath('../src'));
addpath('../data');


load('unit_sphere_4096_noisy_pts.mat');

marker = '+';
linewidth = 2;
    
h = figure;
subplot(121);
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);

scatter3(P(:,1),P(:,2),P(:,3),ones(size(P,1),1),C,marker,'LineWidth',linewidth);
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1],'FontSize',16);
view(3);
title('Noisy point set','Color',[1 1 1],'FontSize',16);

k = 51;
nb_it = 2;
P = smooth_point_set(P,k,nb_it);

subplot(122);
scatter3(P(:,1),P(:,2),P(:,3),ones(size(P,1),1),C,marker,'LineWidth',linewidth);
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1],'FontSize',16);
view(3);
title('Smoothed point set','Color',[1 1 1],'FontSize',16);