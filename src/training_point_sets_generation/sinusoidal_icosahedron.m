function P = sinusoidal_icosahedron(nb_samples, w, random_sampling)
%% sinusoidal_icosahedron : function to compute
% and save a sinusoidal icosahedron point set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
%%% Syntax
%
% sinusoidal_icosahedron;
% sinusoidal_icosahedron(nb_samples);
% sinusoidal_icosahedron(nb_samples, w);
% sinusoidal_icosahedron(nb_samples, w, random_sampling);
% P = sinusoidal_icosahedron(nb_samples, w, random_sampling);
%
%
%%% Description
%
% sinusoidal_icosahedron computes a sinusoidal icosahedron
% point set with parameters nb_samples = 32, w = 1 by default.
%
% sinusoidal_icosahedron(nb_samples) samples at the value nb_samples
% each icosahedron basis triangle (20).
%
% sinusoidal_icosahedron(nb_samples, w) uses the given shape parameter
% w.
%
% sinusoidal_icosahedron(nb_samples, w, random_sampling) generates random samples when random_sampling = true/1.
%
% P = sinusoidal_icosahedron(nb_samples, w)
% stores the point set coordinates in P.
%
%
%%% Input arguments
%
% - nb_samples : positive integer scalar double, nb_samples > 2. Optional.
%   
% - w : real scalar double, the shape parameter.
%       Remarkable value : w = 0 gives a sphere. Optional.
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
        w = 1;
        if nargin < 1
            nb_samples = 32;
        end
    end
end


%% Body
phi_n = 0.5*(1+sqrt(5));

Rmz = [cos(0.4*pi) -sin(0.4*pi) 0;...
       sin(0.4*pi) cos(0.4*pi) 0;...
       0 0 1];

centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));
a = 2/sqrt(phi_n*sqrt(5)); % edge length

% Icosahedron points coordinates
% 1st equilateral triangle
V0 = [0 0 1]';
V1 = [sin(centre_angle) 0 cos(centre_angle)]';
V2 = Rmz*V1;

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

% (12) point set coordinates vector
U0 = Rmz*V2;
U1 = Rmz^2*V2;
U2 = Rmz^3*V2;
U3 = Rmz*V5;
U4 = Rmz^2*V5;
U5 = Rmz^3*V5;

V = [V0 V1 V2 U0 U1 U2 V3 V4 V5 U3 U4 U5];

% Bidirectional (u,v) sampling + compute corresponding squared distances vector
U0 = sample_triangle(V0,V1,V2,nb_samples,random_sampling);
U0 = U0';

% Replicate / rotation -> upper crown
Uu = U0;

for k = 1:4    
    Uu = cat(2,Uu,Rmz^k*U0);    
end

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

U1 = sample_triangle(V3,V4,V5,nb_samples,random_sampling);
U1 = U1';

% Lower crown
Ul = U1;

for k = 1:4
    Ul = cat(2,Ul,Rmz^k*U1);    
end

U = [Uu Ul];

% 1st belt triangle
V6 = V1;
V7 = V2;
V8 = Rmz^2*V5;

U2 = sample_triangle(V6,V8,V7,nb_samples,random_sampling);
U2 = U2';

% 2nd belt triangle
V9  = -V6;
V10 = -V7;
V11 = -V8;

U3 = sample_triangle(V9,V10,V11,nb_samples,random_sampling);
U3 = U3';

% Full belt = centre crown
U4 = [U2 U3];
Uc = U4;

for k = 1:4
    Uc = cat(2,Uc,Rmz^k*U4);    
end

Sinico = [U Uc]; 

% Spherical coordinates
X = Sinico(1,:);
Y = Sinico(2,:);
Z = Sinico(3,:);

X = X(:);
Y = Y(:);
Z = Z(:);
P = [X Y Z];

radius = sqrt(sum(P.^2,2));
P = P./radius;
P = unique(P,'rows');

f = [];
coeff = [];

for i = 1:size(P,1)
    
    % Closest point
    [~,min_dst] = closest_point(P(i,:)',V);        
    
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


end % sinusoidal_icosahedron


%% Closest_point subfunction
function [P0, min_dst] = closest_point(Point, P)


dst_vect = sqrt(sum((P-repmat(Point,[1,size(P,2)])).^2,1));
f = find(dst_vect == min(dst_vect));
min_dst = dst_vect(f(1));
P0 = P(:,f(1));


end % closest_point


%% sample_triangle subfunction
function P = sample_triangle(V1, V2, V3, nb_steps, random_sampling)
%
% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.


% Create sampling grid
Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

P = zeros(sum(1:nb_steps+1),Ndim);

if ~random_sampling
    
    nu = u / norm(u);
    nv = v / norm(v);
    
    stepu = norm(u) / nb_steps;
    stepv = norm(v) / nb_steps;
    
    k = 1;
    
    % Sampling & points generation
    for m = 0:nb_steps
        
        for n = 0:nb_steps
            
            if m+n <= nb_steps % in (V1,V2,V3) triangle conditions ; indices # nb segments
                
                % translation vector
                tv = m*stepu*nu + n*stepv*nv;
                P(k,:) = (V1 + tv)';
                k = k+1;
                
            end
            
        end
        
    end
    
else % if random_sampling
   
    rand_coeff_vect = rand(2,size(P,1)-3);
    f = sum(rand_coeff_vect,1) > 1;
    rand_coeff_vect(:,f) = 1 - rand_coeff_vect(:,f);
    
    % Translation vectors
    TV = repmat(rand_coeff_vect(1,:),[Ndim 1]).*repmat(u,[1 size(rand_coeff_vect,2)]) + ...
         repmat(rand_coeff_vect(2,:),[Ndim 1]).*repmat(v,[1 size(rand_coeff_vect,2)]);        
    
    for k = 1:size(TV,2)
        
        % translation vector
        tv = TV(:,k);                
        P(k,:) = (V1 + tv)';
        
    end
    
    % Add the triangle vertices
    P(end-2,:) = V1';
    P(end-1,:) = V2';
    P(end,:)   = V3';        
    
end


end % sample_triangle