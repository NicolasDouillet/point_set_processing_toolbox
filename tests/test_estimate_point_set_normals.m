% test estimate_point_set_normals

clc;

addpath(genpath('../src'));
addpath('../data');


% % Random unit sphere
% n = 2048;
% V = 2*(rand(n,3) - 0.5);
% V = V ./ sqrt(sum(V.^2,2));


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
N = estimate_point_set_normals(V,nb_ngb,'norm');
plot_point_set(V);
quiver3(V(:,1),V(:,2),V(:,3),N(:,1),N(:,2),N(:,3),'Color',[0 1 1],'LineWidth',1), hold on;
hidden off;
ax = gca;
ax.Clipping= 'off';