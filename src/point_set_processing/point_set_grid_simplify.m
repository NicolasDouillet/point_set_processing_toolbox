function P_out = point_set_grid_simplify(P_in, grid_step, mode)
%% point_set_grid_simplify : function to grid simplify one given point set (P_in).
%
%%% Author : nicolas.douillet (at) free.fr, 2019-2024.
%
%
%%% Input arguments
%
%          [ |  |  |]
% - P_in = [Px Py Pz], real matrix double, the input point set, size(P_in) = [nb_points_in,3].
%          [ |  |  |]
%
% - grid_step : positive integer scalar double the step of the grid.
%
% - mode : character string in the set {*'exact','rounded'}. Case insensitive.
%
%
%%% Output arguments
%
%           [ |  |  |]
% - P_out = [Px Py Pz], real matrix double, the output point set, size(P_out) = [nb_points_out,3].
%           [ |  |  |]


%% Input parsing
assert(nargin > 1, 'Not enough input arguments.');
assert(nargin < 4, 'Too many input arguments.');

if nargin < 3 
    
    mode = 'exact';   
    
else
    
    assert(ischar(mode) && (strcmpi(mode,'exact') || strcmpi(mode,'rounded')),...
           'Input argument mode must be a character string in the set : {''exact'', ''EXACT'',''rounded'',''ROUNDED''}.');
       
end


%% Body
% Zeros padding in 2D case
two_dimensions = false;

if size(P_in,2) == 2
    
    two_dimensions = true;
    P_in = cat(2,P_in,zeros(size(P_in,1),1));
    
end

P_in = unique(P_in,'rows');


% Exact mode
if strcmpi(mode,'exact')

    X_min = min(P_in(:,1));
    X_max = max(P_in(:,1));
    Y_min = min(P_in(:,2));
    Y_max = max(P_in(:,2));
    Z_min = min(P_in(:,3));
    Z_max = max(P_in(:,3));
    
    X_nb_cells = ceil((X_max-X_min)/grid_step);
    Y_nb_cells = ceil((Y_max-Y_min)/grid_step);
    Z_nb_cells = ceil((Z_max-Z_min)/grid_step);
    
    xG_grid = X_min + 0.5*grid_step + grid_step*(0:X_nb_cells-1);
    yG_grid = Y_min + 0.5*grid_step + grid_step*(0:Y_nb_cells-1);
    zG_grid = Z_min + 0.5*grid_step + grid_step*(0:Z_nb_cells-1);
    
    if two_dimensions; zG_grid = Z_min; end
    
    XYZ_cell_idx = ceil((P_in - repmat([X_min, Y_min, Z_min],[size(P_in,1),1]))/grid_step);        
    XYZ_cell_idx(XYZ_cell_idx == 0) = 1;
    
    idx_grid = unique(XYZ_cell_idx,'rows');
    
    P_out = zeros(size(idx_grid,1),3);
    
    for k = 1:size(idx_grid,1)
        
        cell_id = idx_grid(k,:);        
        point_cluster = P_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);                
        G = repmat([xG_grid(1,cell_id(1,1)),yG_grid(1,cell_id(1,2)),zG_grid(1,cell_id(1,3))],[size(point_cluster,1),1]);                
        sq_dst_mat = sum((point_cluster - G).^2,2);
        P_out(k,:) = point_cluster(sq_dst_mat == min(sq_dst_mat),:);
        
    end
    
    
% Rounded mode 
elseif strcmpi(mode,'rounded') % fastest mode
    
    P_out = grid_step*round(P_in/grid_step);
    P_out = unique(P_out,'rows');
    
end


end % point_set_grid_simplify