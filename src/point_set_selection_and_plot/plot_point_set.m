function [] = plot_point_set(P, C)
%% plot_point_set : function to plot the point set (P) in a Matlab (R) figure.
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
%       [| | |]
% - C = [R G B], integer matrix double, the color vector, size(C) = [nb_points,3]. Optional argument.
%       [| | |]
%


%% Input parsing and parameter default values
nb_pts = size(P,1);
marker = '+';
linewidth = 2;

if nargin  < 2
    
    C = repmat([1 1 0],[nb_pts,1]);
    
end
    

%% Body
S = ones(nb_pts,1);

figure;
set(gcf,'Color',[0 0 0]);
scatter3(P(:,1),P(:,2),P(:,3),S,C,marker,'LineWidth',linewidth);
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal;
ax = gca;
set(ax,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
ax.Clipping = 'off';


end % plot_point_set