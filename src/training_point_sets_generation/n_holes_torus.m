function P = n_holes_torus(a, b, Rho, r, n, nb_samples, option_random_sampling)
%% n_holes_torus : function to compute and save a n-holes torus.
%
% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
%%% Input arguments
%
% - a (double)
%
% - b (double)
%
% - Rho (large radius, double)
%
% - r (small radius, double)
%
% - n (nb holes, integer)
%
% - nb_samples (integer)
%
% - option_isotropic (bool)
%
%
%%% Output argument
%
%       [|  |  | ]
% - P = [Px Py Pz], real matrix double, the point set. Size(P) = [nb_points,3].
%       [|  |  | ]


%% Body
X = @(u,v)a*(Rho*cos(v)+r*sin(u).*cos(v));
Y = @(u,v)b*(Rho*sin(v)+r*sin(u).*sin(v));
Z = @(u,v)r*cos(u);

% Sampling parameters
range_u = [0 2*pi floor(0.5*nb_samples)+1];
range_v = [0 2*pi nb_samples];

[M,~,v] = torus_homeo_sfc_isotropic_splg(X,Y,Z,range_u,range_v,option_random_sampling);

% Half torus : sort range_v
f = v < pi;
P = M(f,:);

% Remove junctions points
e = [];

for k = 1:size(P,1)
    
    % ref point
    if k <= floor(0.5*size(P,1))
        
        ref_pt = [Rho 0 0];
        
    elseif k > floor(0.5*size(P,1))
        
        ref_pt = [-Rho 0 0];
        
    end
    
    % vector to data
    vect_2_data = P(k,:) - ref_pt;
    
    % Sorting points to keep
    if abs(vect_2_data(1,3)/vect_2_data(1,2)) <= tan(pi/(n+1))
        
        e = cat(1,e,k);
        
    end
    
end

P = P(e,:);
alpha = pi/n;

Rmx = [1          0           0;
       0 cos(alpha) -sin(alpha);
       0 sin(alpha)  cos(alpha)];

for j = 1:n
    
    R = Rmx^j;
    P = cat(1,P,(R*P')'); 
    
end

P = unique(P,'rows');


end % n_holes_torus