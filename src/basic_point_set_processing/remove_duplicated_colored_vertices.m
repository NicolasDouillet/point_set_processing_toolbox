function [V_out, C_out] = remove_duplicated_colored_vertices(V_in, C_in)
% remove_duplicated_colored_vertices : function to remove
% duplicated colored vertices from the point set.
%
% Author : nicolas.douillet (at) free.fr, 2020-2024.
%
%
% Input arguments
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], real matrix double, the input point set, size(V_in) = [nb_input_vertices,3].
%          [ |    |    |  ]
%
% - C_in : real positive matrix double, the input vertices color set. size(C_in,1) = size(V_in,1).
%
% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], real matrix double, the output point set, size(V_out) = [nb_output_vertices,3],
%           [  |     |     |  ]
%
%           where nb_output_vertices = nb_input_vertices - nb_duplicata.
%
% - C_out : real positive matrix double, the output vertices color set. size(C_out,1) = size(V_out,1).


tol = eps;
[V_out,idx] = uniquetol(V_in,tol,'ByRows',true);
C_out = C_in(idx);


end