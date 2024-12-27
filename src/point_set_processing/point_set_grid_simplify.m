function [P_out, C_out, N_out] = point_set_grid_simplify(P_in, grid_step, mode, C_in, N_in)
%% point_set_grid_simplify : function to grid simplify one given point set (P_in).
%
%%% Author : nicolas.douillet (at) free.fr, 2019-2024.
%
%
%%% Input arguments
%
%          [ |  |  |]
% - P_in = [Px Py Pz], real matrix double, the input point set, size(P_in) = [nb_points_in,3]. Mandatory argument.
%          [ |  |  |]
%
% - grid_step : positive integer scalar double the step of the grid. Mandatory argument.
%
% - mode : character string in the set {*'exact','rounded'}. Case insensitive. Optional argument.
%
%          [ |    |    |  ]
% - C_in = [R_in G_in B_in], integer matrix double, the input color set, size(C_in) = [nb_input_points,3]. Optional argument.
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - N_in = [Nx_in Ny_in Nz_in], real matrix double, the input normals set, size(N_in) = [nb_input_points,3]. Optional argument.
%          [  |     |     |  ]
%
%
%%% Output arguments
%
%           [ |  |  |]
% - P_out = [Px Py Pz], real matrix double, the output point set, size(P_out) = [nb_points_out,3]. Mandatory argument.
%           [ |  |  |]
%
%           [ |     |     |   ]
% - C_out = [R_out G_out B_out], integer matrix double, the output color set, size(C_out) = [nb_output_points,3]. Optional argument.
%           [ |     |     |   ]
%
%           [  |      |      |   ]
% - N_out = [Nx_out Ny_out Nz_out], real matrix double, the output normals set, size(N_in) = [nb_output_points,3]. Optional argument.
%           [  |      |      |   ]


%% Input parsing
assert(nargin > 1, 'Not enough input arguments.');
assert(nargin < 6, 'Too many input arguments.');

if nargin < 3 
    
    mode = 'exact';   
    
else
    
    assert(ischar(mode) && (strcmpi(mode,'exact') || strcmpi(mode,'rounded')),...
           'Input argument mode must be a character string in the set : {''exact'', ''EXACT'',''rounded'',''ROUNDED''}.');
       
end


%% Body
% Zeros padding in 2D case
two_dimensions = false;

if size(P_in,2) < 3
    
    two_dimensions = true;
    P_in = cat(2,P_in,zeros(size(P_in,1),1));
    
end


[P_in,id] = unique(P_in,'rows'); % TODO : with tolerance

if nargin > 3
    
    C_in = C_in(id,:);
    
    if nargin > 4
    
        N_in = N_in(id,:);
        
    end
    
end


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
    
    if nargin < 4
        
        for k = 1:size(idx_grid,1)
            
            cell_id = idx_grid(k,:);
            point_cluster = P_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);
            G = repmat([xG_grid(1,cell_id(1,1)),yG_grid(1,cell_id(1,2)),zG_grid(1,cell_id(1,3))],[size(point_cluster,1),1]);
            sq_dst_mat = sum((point_cluster - G).^2,2);
            P_out(k,:) = point_cluster(sq_dst_mat == min(sq_dst_mat),:);
            
        end
        
        
    elseif nargin > 3
        
        C_out = zeros(size(idx_grid,1),3);
        
        for k = 1:size(idx_grid,1)
            
            cell_id = idx_grid(k,:);
            point_cluster = P_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);
            color_cluster = C_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);
            G = repmat([xG_grid(1,cell_id(1,1)),yG_grid(1,cell_id(1,2)),zG_grid(1,cell_id(1,3))],[size(point_cluster,1),1]);
            sq_dst_mat = sum((point_cluster - G).^2,2);
            P_out(k,:) = point_cluster(sq_dst_mat == min(sq_dst_mat),:);
            C_out(k,:) = color_cluster(sq_dst_mat == min(sq_dst_mat),:);
            
        end
        
        if nargin > 4
            
            N_out = zeros(size(idx_grid,1),3);
            
            for k = 1:size(idx_grid,1)
            
            cell_id = idx_grid(k,:);
            point_cluster = P_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);
            color_cluster = C_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);
            normal_cluster = N_in(all(any(bsxfun(@eq,XYZ_cell_idx,cell_id),3),2),:);
            G = repmat([xG_grid(1,cell_id(1,1)),yG_grid(1,cell_id(1,2)),zG_grid(1,cell_id(1,3))],[size(point_cluster,1),1]);
            sq_dst_mat = sum((point_cluster - G).^2,2);
            P_out(k,:) = point_cluster(sq_dst_mat == min(sq_dst_mat),:);
            C_out(k,:) = color_cluster(sq_dst_mat == min(sq_dst_mat),:);
            N_out(k,:) = normal_cluster(sq_dst_mat == min(sq_dst_mat),:);
            
            end
            
        end
        
    end
    
    
% Rounded mode 
elseif strcmpi(mode,'rounded') % fastest mode
    
    P_out = grid_step*round(P_in/grid_step);
    [P_out,id] = unique(P_out,'rows');
    
    if nargin > 3
        
        C_out = C_in;
        C_out = C_out(id,:);
        
        if nargin > 4
            
            N_out = N_in;
            N_out = N_out(id,:);
            
        end
        
    end
    
end


end % point_set_grid_simplify