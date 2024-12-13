function P_out = remove_points(P_set, P_in, mode)
% remove_points : function to remove points from the point set.
%
% Author & support : nicolas.douillet (at) free.fr, 2020-2024.
%
%
% Input arguments
%
%           [| | |]
% - P_set = [X Y Z], real matrix double, the point set to remove, size(P_set) = [nb_new_points,3],
%           [| | |]
%                    or positive integer row vector double, the index list of the points to remove.
%
%          [ |    |    |  ]
% - P_in = [X_in Y_in Z_in], real matrix double, the input point set, size(P_in) = [nb_input_points,3].
%          [ |    |    |  ]
%
% - mode : character string in the set {'indices','explicit','INDICES','EXPLICIT'}. Explicit mode corresponds to the
%          case where P_set is made of points [X Y Z] coordinates. Case insensitive.
%
%
% Output arguments
%
%           [  |     |     |  ]
% - P_out = [X_out Y_out Z_out], real matrix double, the output point set, size(P_out) = [nb_output_points,3],
%           [  |     |     |  ]
%
%           with nb_output_points = nb_input_points + nb_new_points.


% Body
% tic;
if nargin  < 3    
    mode = 'indices';    
else    
    assert((strcmpi(mode,'indices') && ismember(1,size(P_set))) || (strcmpi(mode,'explicit') && size(P_set,2) == 3),...
           'mode value must be either set to ''indices'' with P_set a one row/column indices vector or to ''explicit'' with size(P_set,2) = 3.');    
end

P_out = P_in;

if strcmpi(mode,'indices') && ismember(1,size(P_set)) || nargin < 4       
    pt_idx_2_rm = P_set;                       
elseif strcmpi(mode,'explicit') && size(P_set,2) == 3       
    pt_idx_2_rm = ismember(P_in,P_set,'rows');    
end

% I Suppress points & triangles
P_out(pt_idx_2_rm,:) = [];
% fprintf('%d points removed in %d seconds.\n',nnz(pt_idx_2_rm),toc);


end % remove_points