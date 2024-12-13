function [P_out, C_out] = remove_duplicated_colored_points(P_in, C_in)
% remove_duplicated_colored_points : function to remove
% duplicated colored points from the point set.
%
%%% Author : nicolas.douillet (at) free.fr, 2020-2024.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
%%% - P_in = [X_in Y_in Z_in], real matrix double, the input point set, size(P_in) = [nb_input_points,3].
%          [ |    |    |  ]
%
% - C_in : real positive matrix double, the input points color set. size(C_in,1) = size(P_in,1).
%
%%% Output arguments
%
%           [  |     |     |  ]
% - P_out = [X_out Y_out Z_out], real matrix double, the output point set, size(P_out) = [nb_output_points,3],
%           [  |     |     |  ]
%
%           where nb_output_points = nb_input_points - nb_duplicata.
%
% - C_out : real positive matrix double, the output points color set. size(C_out,1) = size(P_out,1).


% Body
tol = eps;
[P_out,idx] = uniquetol(P_in,tol,'ByRows',true);
C_out = C_in(idx);


end % remove_duplicated_colored_points