function V = sinusoidal_dodecahedron(sampling, w)
%% sinusoidal_dodecahedron : function to compute
% and save a sinusoidal dodecaahedron point set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
% Syntax
%
% sinusoidal_dodecahedron;
% sinusoidal_dodecahedron(sampling);
% sinusoidal_dodecahedron(sampling, w);
% V = sinusoidal_dodecahedron(sampling, w);
%
%
%%% Description
%
% sinusoidal_dodecahedron computes a sinusoidal icosahedron point set
% with parameters sampling = 60, w = 1 by default.
%
% sinusoidal_dodecahedron(sampling) samples at the value sampling
% each icosahedron basis triangle (20).
%
% sinusoidal_dodecahedron(sampling, w) uses the given shape parameter w.
%
% [V, T] = sinusoidal_dodecahedron(sampling, w) stores
% the point set coordinates in V.
%
%
%%% Input arguments
%
% - sampling : positive integer scalar double, sampling > 2.
%              Remarkable value : sampling = 3 gives an icosahedron. Optional.
%   
% - w : real scalar double, the shape parameter.
%       Remarkable value : w = 0 gives a geoid. Optional.
%
%
%%% Output arguments
%
%       [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the point set. Size(V) = [nb_points,3].
%       [|  |  | ]
%
%
%%% Example #1 : default parameters values
% sinusoidal_dodecahedron;
%
%
%%% Example #2 : minimum sampling step
% sinusoidal_dodecahedron(6);
%
%
%%% Example #3 : negative shape parameter value
% sinusoidal_dodecahedron(60,-1);
%
%
%%% Example #4 : large shape parameter value
% sinusoidal_dodecahedron(60,3);


%% Input parsing and default values
if nargin < 2
    w = 1;
    if nargin < 1
        sampling = 60;
        
    end
end


%% Body
phi_n = 0.5*(1+sqrt(5));

Mrz = [cos(0.4*pi) -sin(0.4*pi) 0;...
       sin(0.4*pi) cos(0.4*pi) 0;...
       0 0 1];

centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));
a = 4/3/sqrt(phi_n*sqrt(5)); % edge length

% Icosahedron points coordinates
% 1st equilateral triangle
V0 = [0 0 1]';
V1 = [sin(centre_angle) 0 cos(centre_angle)]';
V2 = Mrz*V1;

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

% (12) points set coordinates vector
U0 = Mrz*V2;
U1 = Mrz^2*V2;
U2 = Mrz^3*V2;
U3 = Mrz*V5;
U4 = Mrz^2*V5;
U5 = Mrz^3*V5;

% Icosahedron face centres
C0 = mean(cat(2,V0,V1,V2),2)';
C1 = mean(cat(2,V0,V2,U0),2)';
C2 = mean(cat(2,V0,U0,U1),2)';
C3 = mean(cat(2,V0,U1,U2),2)';
C4 = mean(cat(2,V0,U2,V1),2)';

C5 = mean(cat(2,V3,V4,V5),2)';
C6 = mean(cat(2,V3,V5,U3),2)';
C7 = mean(cat(2,V3,U3,U4),2)';
C8 = mean(cat(2,V3,U4,U5),2)';
C9 = mean(cat(2,V3,U5,V4),2)';

C10 = mean(cat(2,V1,V2,U4),2)';
C11 = mean(cat(2,V2,U4,U5),2)';
C12 = mean(cat(2,V2,U0,U5),2)';
C13 = mean(cat(2,U0,U5,V4),2)';
C14 = mean(cat(2,U0,U1,V4),2)';

C15 = mean(cat(2,U1,V4,V5),2)';
C16 = mean(cat(2,U1,U2,V5),2)';
C17 = mean(cat(2,U2,V5,U3),2)';
C18 = mean(cat(2,U2,V1,U3),2)';
C19 = mean(cat(2,V1,U3,U4),2)';

C = cat(1,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19)';

% Bidirectional (u,v) sampling + compute corresponding squared distances vector
U0 = sample_triangle(V0,V1,V2,floor(sampling/3));
U0 = U0';

% Replicate / rotation -> upper crown
Uu = U0;

for k = 1:4    
    Uu = cat(2,Uu, Mrz^k*U0);    
end

U1 = sample_triangle(V3,V4,V5,floor(sampling/3));
U1 = U1';

% Lower crown
Ul = U1;

for k = 1:4
    Ul = cat(2,Ul,Mrz^k*U1);    
end

U = [Uu Ul];

% 1st belt triangle
V6 = V1;
V7 = V2;
V8 = Mrz^2*V5;

U2 = sample_triangle(V6,V8,V7,floor(sampling/3));
U2 = U2';

% 2nd belt triangle
V9  = -V6;
V10 = -V7;
V11 = -V8;

U3 = sample_triangle(V9,V10,V11,floor(sampling/3));
U3 = U3';

% Full belt = centre crown
U4 = [U2 U3];
Uc = U4;

for k = 1:4
    Uc = cat(2,Uc,Mrz^k*U4);    
end

Sinico = [U Uc]; 

% Spherical coordinates
X = Sinico(1,:);
Y = Sinico(2,:);
Z = Sinico(3,:);

X = X(:);
Y = Y(:);
Z = Z(:);
V = [X Y Z];
V = unique(V,'rows');

for i = 1:size(V,1)
    
    radius = sqrt(sum(V(i,:).^2));
    
    if radius ~= 1
        
        V(i,:) = V(i,:) / radius;
        
    end
    
end

f = [];
coeff = [];

for i = 1:size(V,1)
    
    % Closest point
    [~, min_dst] = closest_point(V(i,:)',C);        
    
    if min_dst < 0.5*a
        
        f = cat(1,f,i);
        coeff = cat(1,coeff,2*pi*min_dst/a);
        
    end
    
end

N = V(f,:);

X_s = N(:,1);
Y_s = N(:,2);
Z_s = N(:,3);

radius = sqrt(sum(X_s.^2+Y_s.^2+Z_s.^2, 2));
theta_c = acos(Z_s ./ radius);
phi_c = atan2(Y_s, X_s);
Rho_s = w*0.5*(1+cos(coeff));

% Arcos option
X_s = X_s + Rho_s.*sin(theta_c).*cos(phi_c);
Y_s = Y_s + Rho_s.*sin(theta_c).*sin(phi_c);
Z_s = Z_s + Rho_s.*cos(theta_c);

N = [X_s Y_s Z_s];
V(f,:) = N;
V = unique(V,'rows');


end % sinusoidal_dodecahedron


%% Closest_point subfunction
function [V0, min_dst] = closest_point(Point, V)


dst_vect = sqrt(sum((V-repmat(Point, [1,size(V,2)])).^2,1));
f = find(dst_vect == min(dst_vect));
min_dst = dst_vect(f(1));
V0 = V(:,f(1));


end % closest_point


%% sample_triangle subfunction
function V = sample_triangle(V1, V2, V3, nbstep)
%
% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.


% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1),Ndim);

nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nbstep;
stepv = norm(v) / nbstep;
k = 1;

% Sampling & points generation
for m = 0:nbstep
    
    for n = 0:nbstep
        
        if m+n <= nbstep % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end


end % sample_triangle