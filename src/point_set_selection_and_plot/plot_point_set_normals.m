function [] = plot_point_set_normals(P, N, C)
%% plot_point_set_normals : function to plot the point set (P) and its normals set in a Matlab (R) figure.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2020-2025.
%
%
%%% Input arguments
%
%       [| | |]
% - P = [X Y Z], real matrix double, the point set, size(P) = [nb_points,3]. Mandatory argument.
%       [| | |]
%
%       [ |  |  |]
% - N = [Nx Ny Nz], real matrix double, the normals set, size(N) = [nb_input_points,3]. Mandatory argument.
%       [ |  |  |]
%
%       [| | |]
% - C = [R G B], integer matrix double, the color vector, size(C) = [nb_points,3]. Optional argument.
%       [| | |]


%% Input parsing and parameter default values
nb_pts = size(P,1);
marker = '+';
linewidth = 2;

if nargin  < 3
    
    C = repmat([1 1 0],[nb_pts,1]);
    
end


%% Body
S = ones(nb_pts,1);

figure;
set(gcf,'Color',[0 0 0]);
scatter3(P(:,1),P(:,2),P(:,3),S,C,marker,'LineWidth',linewidth), hold on;
quiver3(P(:,1),P(:,2),P(:,3),N(:,1),N(:,2),N(:,3),'Color',[0 1 1],'LineWidth',1), hold on;
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
ax = gca;
set(ax,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
ax.Clipping = 'off';


end % plot_point_set_normals