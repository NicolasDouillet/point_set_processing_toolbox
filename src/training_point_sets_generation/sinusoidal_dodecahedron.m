function P = sinusoidal_dodecahedron(nb_samples, w, random_sampling)
%% sinusoidal_dodecahedron : function to compute
% and save a sinusoidal dodecaahedron point set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
%%% Syntax
%
% sinusoidal_dodecahedron;
% sinusoidal_dodecahedron(nb_samples);
% sinusoidal_dodecahedron(nb_samples, w);
% sinusoidal_dodecahedron(nb_samples, w, random_sampling);
% P = sinusoidal_dodecahedron(nb_samples, w, random_sampling);
%
%
%%% Description
%
% sinusoidal_dodecahedron computes a sinusoidal icosahedron point set
% with parameters nb_samples = 180, w = 3 by default.
%
% sinusoidal_dodecahedron(nb_samples) samples at the value nb_samples
% each icosahedron basis triangle (20).
%
% sinusoidal_dodecahedron(nb_samples, w) uses the given shape parameter w.
%
% % sinusoidal_dodecahedron(nb_samples, w, random_sampling) generates random samples if random_sampling = true/1.
%
% P = sinusoidal_dodecahedron(nb_samples, w, random_sampling) stores the point set coordinates in P.
%
%
%%% Input arguments
%
% - nb_samples : positive integer scalar double, nb_samples > 2. Optional.
%   
% - w : real scalar double, the shape parameter.
%       Remarkable value : w = 0 gives a geoid. Optional.
%
% - random_sampling : logical false*/true | 0*/1. Optional.
%
%
%%% Output arguments
%
%       [|  |  | ]
% - P = [Px Py Pz], real matrix double, the point set. Size(P) = [nb_points,3].
%       [|  |  | ]


%% Input parsing and default values
if nargin < 3
    random_sampling = false;
    if nargin < 2
        w = 3;
        if nargin < 1
            nb_samples = 180;
        end
    end
end


%% Body
phi_n = 0.5*(1+sqrt(5));

Rmz = [cos(0.4*pi) -sin(0.4*pi) 0;...
       sin(0.4*pi) cos(0.4*pi) 0;...
       0 0 1];

centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));
a = 4/3/sqrt(phi_n*sqrt(5)); % edge length

% Icosahedron points coordinates
% 1st equilateral triangle
V0 = [0 0 1]';
V1 = [sin(centre_angle) 0 cos(centre_angle)]';
V2 = Rmz*V1;

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

% (12) points set coordinates vector
U0 = Rmz*V2;
U1 = Rmz^2*V2;
U2 = Rmz^3*V2;
U3 = Rmz*V5;
U4 = Rmz^2*V5;
U5 = Rmz^3*V5;

% Icosahedron face centres
C0 = mean([V0,V1,V2],2)';
C1 = mean([V0,V2,U0],2)';
C2 = mean([V0,U0,U1],2)';
C3 = mean([V0,U1,U2],2)';
C4 = mean([V0,U2,V1],2)';

C5 = mean([V3,V4,V5],2)';
C6 = mean([V3,V5,U3],2)';
C7 = mean([V3,U3,U4],2)';
C8 = mean([V3,U4,U5],2)';
C9 = mean([V3,U5,V4],2)';

C10 = mean([V1,V2,U4],2)';
C11 = mean([V2,U4,U5],2)';
C12 = mean([V2,U0,U5],2)';
C13 = mean([U0,U5,V4],2)';
C14 = mean([U0,U1,V4],2)';

C15 = mean([U1,V4,V5],2)';
C16 = mean([U1,U2,V5],2)';
C17 = mean([U2,V5,U3],2)';
C18 = mean([U2,V1,U3],2)';
C19 = mean([V1,U3,U4],2)';

C = cat(1,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19)';


if random_sampling
    
    P = random_unit_sphere_point_set(0.5*nb_samples^2);
    P = P./sqrt(sum(P.^2,2));
    
else
    
    X = @(u,v)sin(u).*cos(v);
    Y = @(u,v)sin(u).*sin(v);
    Z = @(u,v)cos(u);
    
    range_u = [0 pi nb_samples/2+1];
    range_v = [0 2*pi nb_samples];
    
    P = sphere_homeo_sfc_isotropic_splg(X,Y,Z,range_u,range_v,false);
    
end

f = [];
coeff = [];

for i = 1:size(P,1)
    
    % Closest point
    [~,min_dst] = closest_point(P(i,:)',C);        
    
    if min_dst < 0.5*a
        
        f = cat(1,f,i);
        coeff = cat(1,coeff,2*pi*min_dst/a);
        
    end
    
end

N = P(f,:);

X_s = N(:,1);
Y_s = N(:,2);
Z_s = N(:,3);

radius = sqrt(sum(X_s.^2+Y_s.^2+Z_s.^2,2));
theta_c = acos(Z_s ./ radius);
phi_c = atan2(Y_s,X_s);
Rho_s = w*0.5*(1+cos(coeff));

% Arcos option
X_s = X_s + Rho_s.*sin(theta_c).*cos(phi_c);
Y_s = Y_s + Rho_s.*sin(theta_c).*sin(phi_c);
Z_s = Z_s + Rho_s.*cos(theta_c);

N = [X_s Y_s Z_s];
P(f,:) = N;
P = unique(P,'rows');


end % sinusoidal_dodecahedron


%% Closest_point subfunction
function [P0, min_dst] = closest_point(Point, P)


dst_vect = sqrt(sum((P-repmat(Point,[1,size(P,2)])).^2,1));
f = find(dst_vect == min(dst_vect));
min_dst = dst_vect(f(1));
P0 = P(:,f(1));


end % closest_point