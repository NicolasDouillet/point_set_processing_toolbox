% test estimate_point_set_normals

clc;

addpath(genpath('../src'));
addpath('../data');


% % Random unit sphere
% n = 2048;
% P = 2*(rand(n,3) - 0.5);
% P = P ./ sqrt(sum(P.^2,2));


filenames = {'torus';...             
             'bowser';...
             'vase';...
             'kitten';...
             'virus_cell_like_surface_light';...
             'sinico';...
             'bowser_MD';...
             'bunny_16k'
             };

fid = 4;
filename = strcat(cell2mat(filenames(fid,1)),'.mat');         
load(filename);


nb_ngb = 4; % >= 4 and adapt function of point set local density 
N = estimate_point_set_normals(P,nb_ngb); % normalized and oriented by default
plot_point_set(P);
quiver3(P(:,1),P(:,2),P(:,3),N(:,1),N(:,2),N(:,3),'Color',[0 1 1],'LineWidth',1), hold on;
hidden off;
ax = gca;
ax.Clipping= 'off';