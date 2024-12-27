function P = create_random_unit_sphere_point_set(n)
%% create_random_unit_sphere_point_set : function to create random noisy sphere in spherical coordinates
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input argument
%
% - n : positive integer scalar double, the number of points / samples in
%       the point set.
%
%
%%% Output argument
%
%       [| | |]
% - P = [X Y Z], real matrix double, the point set, size(P) = [nb_points,3].
%       [| | |]


%% Body
M = rand(3,n);

M(1,:) = 1 + 0.2*(M(1,:) - 0.5);
M(2,:) = pi*M(2,:);
M(3,:) = 2*pi*M(3,:);

rho   = M(1,:);
theta = M(2,:);
phi   = M(3,:);

X = rho.*sin(theta).*cos(phi);
Y = rho.*sin(theta).*sin(phi);
Z = rho.*cos(theta);

P = cat(2,X',Y',Z');


end % create_random_unit_sphere_point_set