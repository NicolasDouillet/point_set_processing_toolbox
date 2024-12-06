function V_out = remove_duplicated_vertices(V_in)
% remove_duplicated_vertices : function to remove
% duplicated vertices from the point set.
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
% Output arguments
%
%           [  |     |     |  ]
% - V_out = [X_out Y_out Z_out], real matrix double, the output point set, size(V_out) = [nb_output_vertices,3],
%           [  |     |     |  ]
%
%           where nb_output_vertices = nb_input_vertices - nb_duplicata.


% Body
% tic;
tol = eps;
V_out = uniquetol(V_in,tol,'ByRows',true);

% fprintf('%d duplicated vertices removed in %d seconds.\n',size(V_in,1)-size(V_out,1),toc);


end