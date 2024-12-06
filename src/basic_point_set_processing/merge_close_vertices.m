function [V_out, T_out] = merge_close_vertices(V_in, T_in, tol)
% merge_close_vertices : function to merge vertices of the point closer
% than the tolerance tol in terms of euclidian distance.
%
% Author & support : nicolas.douillet (at) free.fr, 2020.
%
%
% Input arguments
%
%          [ |    |    |  ]
% - V_in = [X_in Y_in Z_in], double matrix, the input point set, size(V_in) = [nb_input_vertices,3].
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix, the input triangulation, size(T_in) = [nb_input_triangles,3].
%          [  |     |     |  ]
%
% - tol : double scalar, the 3D distance tolerance.
%
%
% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], double matrix, the output point set, size(V_out) = [nb_output_vertices,3].
%           [  |     |     |  ]
%
%           [  |      |      |   ]
% - T_out = [i1_out i2_out i3_out], positive integer matrix, the output triangulation, size(T_out) = [nb_output_triangles,3].
%           [  |      |      |   ]


% tic;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);
% fprintf('merge_close_vertices request executed in %d seconds.\n',toc);


end % merge_close_vertices