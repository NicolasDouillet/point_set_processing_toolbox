function V_out = point_set_grid_simplify(V_in, grid_step, version)
%% point_set_grid_simplify : function to grid simplify one given point set (V_in).
%
%%% Author : nicolas.douillet (at) free.fr, 2019-2024.
%
%
%%% Syntax
%
% V_out = point_set_grid_simplify(V_in, grid_step);
%
% V_out = point_set_grid_simplify(V_in, grid_step, version);
%
%
%%% Description
%
% V_out = point_set_grid_simplify(V_in, grid_step) simplifies the input
% point set V_in through the 3D cubic grid of step grid_step and return
% the result V_out.
%
%%% V_out = point_set_grid_simplify(V_in, grid_step, version) follows version option
% *'exact'/'rounded' for the simplification process. In the exact version,
% every vertex of V_out also belong to V_in set : for each cell of the grid
% the closest to the cell centre is kept, whereas rounded version is an
% approximation (for each cell of the grid, the corresponding vertex cluster is
% replaced with the cell centre) though sometimes usefull because very fast.
%
%
%%% See also COLON, LINSPACE
%
%
%%% Input arguments
%
%          [ |  |  |]
% - V_in = [Vx Vy Vz], real matrix double, the input point set, size(V_in) = [nb_vertices_in,3].
%          [ |  |  |]
%
% - grid_step : positive integer scalar double the step of the grid.
%
% - version : character string in the set {*'exact','rounded'}. Case insensitive.
%
%
%%% Output arguments
%
%           [ |  |  |]
% - V_out = [Vx Vy Vz], real matrix double, the output point set, size(V_out) = [nb_vertices_out,3].
%           [ |  |  |]
%
%
%%% Example : unit ball sampling
% N = 1e4;
% X = 2*(rand(N,1)-0.5);
% Y = 2*(rand(N,1)-0.5);
% Z = 2*(rand(N,1)-0.5);
% 
% Rho = X.^2 + Y.^2 + Z.^2;
% i = Rho <= 1;
% X = X(i);
% Y = Y(i);
% Z = Z(i);
% 
% M = [X, Y, Z];
% grid_step = 0.2;
% N = point_set_grid_simplify(M, grid_step);
%
% figure;
% subplot(121);
% plot3(X,Y,Z,'.','Color',[0 0 1],'Linewidth',2), hold on;
% axis equal, view(2);
% title({'Input raw point set',cat(2,num2str(size(M,1)),' vertices')});
%
% subplot(122);
% plot3(N(:,1), N(:,2), N(:,3),'.','Color',[0 0 1],'Linewidth',2), hold on;
% axis equal, view(2);
% title({'Output grid sampled point set',cat(2,'grid step = ', num2str(grid_step),' ; ',num2str(size(N,1)),' vertices')});


%% Input parsing
assert(nargin > 1, 'Not enough input arguments.');
assert(nargin < 4, 'Too many input arguments.');

if nargin < 3 
    
    version = 'exact';   
    
else
    
    assert(ischar(version) && (strcmpi(version,'exact') || strcmpi(version,'rounded')),...
           'Input argument version must be a character string in the set : {''exact'', ''EXACT'',''rounded'',''ROUNDED''}.');
       
end

% Input data format
assert(isnumeric(grid_step) && isreal(grid_step) && grid_step > 0, 'Input argument grid_step must be a positive integer.');
assert(isnumeric(V_in) && isreal(V_in) && size(V_in,1) > 0 && size(V_in,2) > 1 && size(V_in,2) < 4, 'Input argument V_in must be an N x 3 or N x 2 real numeric matrix (N : vertex number).')


%% Body
% Zeros padding in 2D case
two_dimensions = false;

if size(V_in,2) == 2
    
    two_dimensions = true;
    V_in = cat(2,V_in,zeros(size(V_in,1),1));
    
end

V_in = unique(V_in,'rows');


% Exact version
if strcmpi(version,'exact')

    X_min = min(V_in(:,1));
    X_max = max(V_in(:,1));
    Y_min = min(V_in(:,2));
    Y_max = max(V_in(:,2));
    Z_min = min(V_in(:,3));
    Z_max = max(V_in(:,3));
    
    X_nb_cells = ceil((X_max-X_min)/grid_step);
    Y_nb_cells = ceil((Y_max-Y_min)/grid_step);
    Z_nb_cells = ceil((Z_max-Z_min)/grid_step);
    
    xG_grid = X_min + 0.5*grid_step + grid_step*(0:X_nb_cells-1);
    yG_grid = Y_min + 0.5*grid_step + grid_step*(0:Y_nb_cells-1);
    zG_grid = Z_min + 0.5*grid_step + grid_step*(0:Z_nb_cells-1);
    
    if two_dimensions; zG_grid = Z_min; end;
    
    XYZ_cell_idx = ceil((V_in - repmat([X_min, Y_min, Z_min],[size(V_in,1),1]))/grid_step);        
    XYZ_cell_idx(XYZ_cell_idx == 0) = 1;
    
    idx_grid = unique(XYZ_cell_idx,'rows');
    
    V_out = zeros(size(idx_grid,1),3);
    
    for k = 1:size(idx_grid,1)
        
        cell_id = idx_grid(k,:);        
        vertex_cluster = V_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);                
        G = repmat([xG_grid(1,cell_id(1,1)),yG_grid(1,cell_id(1,2)),zG_grid(1,cell_id(1,3))],[size(vertex_cluster,1),1]);                
        sq_dst_mat = sum((vertex_cluster - G).^2,2);
        V_out(k,:) = vertex_cluster(sq_dst_mat == min(sq_dst_mat),:);
        
    end
    
    
% Rounded version 
elseif strcmpi(version,'rounded') % fastest version
    
    V_out = grid_step*round(V_in/grid_step);
    V_out = unique(V_out,'rows');
    
end


end % point_set_grid_simplify