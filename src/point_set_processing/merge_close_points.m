function P_out = merge_close_points(P_in, tol)
% merge_close_points : function to merge points of the point closer
% than the tolerance tol in terms of euclidian distance.
%
%%% Author & support : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
% - P_in = [X_in Y_in Z_in], double matrix, the input point set, size(P_in) = [nb_input_points,3].
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix, the input triangulation, size(T_in) = [nb_input_triangles,3].
%          [  |     |     |  ]
%
% - tol : double scalar, the 3D distance tolerance.
%
%
%%% Output arguments
%
%           [  |     |     |  ]
% - P_out = [X_out Y_out Z_out], double matrix, the output point set, size(P_out) = [nb_output_points,3].
%           [  |     |     |  ]
%
%           [  |      |      |   ]
% - T_out = [i1_out i2_out i3_out], positive integer matrix, the output triangulation, size(T_out) = [nb_output_triangles,3].
%           [  |      |      |   ]


% Body
% tic;
P_out = uniquetol(P_in,tol,'ByRows',true);
% fprintf('merge_close_points request executed in %d seconds.\n',toc);


end % merge_close_points