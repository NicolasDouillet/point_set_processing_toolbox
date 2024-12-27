function [P_out, C_out, N_out] = merge_close_points(P_in, tol, C_in, N_in)
%% merge_close_points : function to merge together points closer
% than the tolerance tol in terms of euclidian distance.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
% - P_in = [X_in Y_in Z_in], double matrix, the input point set, size(P_in) = [nb_input_points,3]. Mandatory argument.
%          [ |    |    |  ]
%
% - tol : double scalar, the 3D distance tolerance. Mandatory argument.
%
%          [ |    |    |  ]
% - C_in = [R_in G_in B_in], integer matrix double, the input color set, size(C_in) = [nb_input_points,3]. Optional argument.
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - N_in = [Nx_in Ny_in Nz_in], real matrix double, the input normals set, size(N_in) = [nb_input_points,3]. Optional argument.
%          [  |     |     |  ]
%
%%% Output arguments
%
%           [  |     |     |  ]
% - P_out = [X_out Y_out Z_out], double matrix, the output point set, size(P_out) = [nb_output_points,3]. Mandatory argument.
%           [  |     |     |  ]
%
%           [ |     |     |   ]
% - C_out = [R_out G_out B_out], integer matrix double, the output color set, size(C_out) = [nb_output_points,3]. Optional argument.
%           [ |     |     |   ]
%
%           [  |      |      |   ]
% - N_out = [Nx_out Ny_out Nz_out], real matrix double, the output normals set, size(N_in) = [nb_output_points,3]. Optional argument.
%           [  |      |      |   ]


%% Body
% tic;
[P_out,id] = uniquetol(P_in,tol,'ByRows',true);

if nargin > 2
   
    C_out = C_in(id,:);
    
    if nargin > 3
        
        N_out = N_in(id,:);
        
    end
    
end

% fprintf('merge_close_points request executed in %d seconds.\n',toc);


end % merge_close_points