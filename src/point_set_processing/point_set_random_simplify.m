function [P_out, C_out, N_out] = point_set_random_simplify(P_in, coeff, C_in, N_in)
%% point_set_random_simplify : function to random simplify one given point set (P_in).
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input arguments
%
%          [ |  |  |]
% - P_in = [Px Py Pz], real matrix double, the input point set, size(P_in) = [nb_points_in,3]. Mandatory argument.
%          [ |  |  |]
%
% - coeff : positive real scalar double, the simplification coefficient. 0 <= coeff <= 1.
%           coeff = x means output P_out will contain 100*x % the number of points in P_in.
%           Default value is coeff = 0.5. Optional argument.
%
%          [ |    |    |  ]
% - C_in = [R_in G_in B_in], integer matrix double, the input color set, size(C_in) = [nb_input_points,3]. Optional argument.
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - N_in = [Nx_in Ny_in Nz_in], real matrix double, the input normals set, size(N_in) = [nb_input_points,3]. Optional argument.
%          [  |     |     |  ]
%
%
%%% Output arguments
%
%           [ |  |  |]
% - P_out = [Px Py Pz], real matrix double, the output point set, size(P_out) = [nb_points_out,3]. Mandatory argument.
%           [ |  |  |]
%
%           [ |     |     |   ]
% - C_out = [R_out G_out B_out], integer matrix double, the output color set, size(C_out) = [nb_output_points,3]. Optional argument.
%           [ |     |     |   ]
%
%           [  |      |      |   ]
% - N_out = [Nx_out Ny_out Nz_out], real matrix double, the output normals set, size(N_in) = [nb_output_points,3]. Optional argument.
%           [  |      |      |   ]


%% Input parsing
if nargin > 1
    
    if ~isnumeric(coeff) | ~isreal(coeff) | coeff < 0 | coeff > 1
        
        error('Simplification coefficient coeff must be a real positive number in the range [0; 1].');
        
    end
    
else % if nargin < 2
    
    coeff = 0.5;
    
end


%% Body
nb_pts = size(P_in,1);
id = randi(nb_pts,1,round(coeff*nb_pts)); 
P_out = P_in(id,:);

if nargin > 2
    
    C_out = C_in(id,:);
    
    if nargin > 3
        
        N_out = N_in(id,:);
        
    end
    
end


end % point_set_random_simplify