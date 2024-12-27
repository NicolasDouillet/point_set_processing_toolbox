function [P_out, C_out, N_out] = remove_duplicated_points(P_in, C_in, N_in)
%% remove_duplicated_points : function to remove
% duplicated points from the point set.
%
%%% Author : nicolas.douillet (at) free.fr, 2020-2024.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
% - P_in = [X_in Y_in Z_in], real matrix double, the input point set, size(P_in) = [nb_input_points,3]. Mandatory argument.
%          [ |    |    |  ]
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
% - P_out = [X_out Y_out Z_out], real matrix double, the output point set, size(P_out) = [nb_output_points,3],
%           [  |     |     |  ]
%
%           where nb_output_points = nb_input_points - nb_duplicata. Mandatory argument.
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
tol = eps;
[P_out,id] = uniquetol(P_in,tol,'ByRows',true);

if nargin > 1
   
    C_out = C_in(id,:);
    
    if nargin > 2
        
        N_out = N_in(id,:);
        
    end
    
end

% fprintf('%d duplicated points removed in %d seconds.\n',size(P_in,1)-size(P_out,1),toc);


end % remove_duplicated_points