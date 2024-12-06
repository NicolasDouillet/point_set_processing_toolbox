function V_out = remove_vertices(V_set, V_in, mode)
% remove_vertices : function to remove vertices from the vertex set.
%
% Author & support : nicolas.douillet (at) free.fr, 2020-2024.
%
%
% Input arguments
%
%           [| | |]
% - V_set = [X Y Z], real matrix double, the vertex set to remove, size(V_set) = [nb_new_vertices,3],
%           [| | |]
%                    or positive integer row vector double, the index list of the vertices to remove.
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], real matrix double, the input point set, size(V_in) = [nb_input_vertices,3].
%          [ |    |    |  ]
%
% - mode : character string in the set {'indices','explicit','INDICES','EXPLICIT'}. Explicit mode corresponds to the
%          case where V_set is made of vertices [X Y Z] coordinates. Case insensitive.
%
%
% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], real matrix double, the output point set, size(V_out) = [nb_output_vertices,3],
%           [  |     |     |  ]
%
%           with nb_output_vertices = nb_input_vertices + nb_new_vertices.


% Body
% tic;
if nargin  < 3    
    mode = 'indices';    
else    
    assert((strcmpi(mode,'indices') && ismember(1,size(V_set))) || (strcmpi(mode,'explicit') && size(V_set,2) == 3),...
           'mode value must be either set to ''indices'' with V_set a one row/column indices vector or to ''explicit'' with size(V_set,2) = 3.');    
end

V_out = V_in;

if strcmpi(mode,'indices') && ismember(1,size(V_set)) || nargin < 4       
    vtx_idx_2_rm = V_set;                       
elseif strcmpi(mode,'explicit') && size(V_set,2) == 3       
    vtx_idx_2_rm = ismember(V_in,V_set,'rows');    
end

% I Suppress vertices & triangles
V_out(vtx_idx_2_rm,:) = [];
% fprintf('%d vertices removed in %d seconds.\n',nnz(vtx_idx_2_rm),toc);


end % remove_vertices