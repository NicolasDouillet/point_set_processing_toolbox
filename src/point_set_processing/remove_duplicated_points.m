function P_out = remove_duplicated_points(P_in)
% remove_duplicated_points : function to remove
% duplicated points from the point set.
%
%%% Author : nicolas.douillet (at) free.fr, 2020-2024.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
% - P_in = [X_in Y_in Z_in], real matrix double, the input point set, size(P_in) = [nb_input_points,3].
%          [ |    |    |  ]
%
%%% Output arguments
%
%           [  |     |     |  ]
% - P_out = [X_out Y_out Z_out], real matrix double, the output point set, size(P_out) = [nb_output_points,3],
%           [  |     |     |  ]
%
%           where nb_output_points = nb_input_points - nb_duplicata.


% Body
% tic;
tol = eps;
P_out = uniquetol(P_in,tol,'ByRows',true);

% fprintf('%d duplicated points removed in %d seconds.\n',size(P_in,1)-size(P_out,1),toc);


end % remove_duplicated_points