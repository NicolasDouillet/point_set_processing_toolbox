function P = Reuleaux_tetrahedron(nb_samples, random_sampling)
%% Reuleaux_tetrahedron : function to generate a Reuleaux tetrahedron point
% set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2017-2025.
%
%
%%% Input arguments
%
% - nb_samples : positive integer scalar double, nb_samples > 2. Default value is nb_samples = 32. Optional.
%
% - random_sampling : logical false*/true | 0*/1. Optional.
%
%
%%% Ouput argument
%
%       [|  |  | ]
% - P = [Px Py Pz], real matrix double, the point set. Size(P) = [nb_points,3].
%       [|  |  | ]


%% Input parsing
if nargin < 2
    
   random_sampling = false ;
   
   if nargin < 1
      
       nb_samples = 32;
       
   end
   
end


%% Body
% I Summits of original tetrahedron (living in the unit sphere R(O,1))
P1 = [0 0 1]';
P2 = [2*sqrt(2)/3 0 -1/3]';
P3 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
P4 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';


edge_length = norm(P1-P2); %  = 2*sqrt(6)/3

P123 = sample_and_curve_triangle(P1,P2,P3,nb_samples,0.85,random_sampling);
P123 = inflate_triangle_sample_from_opposite_vertex(P123,P4',edge_length);

% Tetrahedron faces rotations
Rmy = @(theta) [cos(theta) 0 -sin(theta);
                0          1  0;
                sin(theta) 0  cos(theta)];

Rmz = @(theta) [cos(theta) -sin(theta) 0;
                sin(theta)  cos(theta) 0;
                0           0          1];
                        
P134 = (Rmz(2*pi/3)*P123')';
P142 = (Rmz(2*pi/3)*P134')';
P234 = (Rmy(-acos(-1/3))*Rmz(pi/3)*P142')';
            
P = [P123; P134; P142; P234];
P = unique(P, 'rows');


end % Reuleaux_tetrahedron


%% sample_and_curve_triangle subfunction
function I = sample_and_curve_triangle(P1, P2, P3, nb_samples, warp, random_sampling)
%
% Authors : nicolas.douillet (at) free.fr, 2017-2024.
%           Gerd Wachsmuth,                     2021.


% Create sampling grid
Ndim = size(P1,1);

% (P0P1, P0P2) base
u = (P2 - P1);
v = (P3 - P1);

I = zeros((nb_samples+1)*(nb_samples+2)/2, Ndim);
    
if ~random_sampling
    
    k = 1;
    
    % Sampling & vertices generation
    for m = 0:nb_samples
        
        for n = 0:nb_samples
            
            if (m+n) <= nb_samples % in (P0,P1,P2) triangle conditions ; indices # nb segments
                
                % Barycentric coordinates.
                l1 = m/nb_samples;
                l2 = n/nb_samples;
                l3 = (nb_samples - n - m)/nb_samples;
                
                % Transform the barycentric coordinates.
                b1 = l1^warp;
                b2 = l2^warp;
                b3 = l3^warp;
                
                % Assure that they still sum up to 1.
                db = (b1 + b2 + b3) - 1;
                b1 = b1 - db*l1;
                b2 = b2 - db*l2;
                % b3 = b3 - db*l3;
                
                % translation vector
                tv = b1*u + b2*v;
                I(k,:) = (P1 + tv)';
                k = k+1;
                
            end
            
        end
        
    end
    
else % if random_sampling
    
    nb_points = sum(1:nb_samples);
    I = zeros(nb_points,Ndim);    
    rand_coeff_vect = rand(2,nb_points-3);
    f = sum(rand_coeff_vect,1) > 1;
    rand_coeff_vect(:,f) = 1 - rand_coeff_vect(:,f);
    
    % Translation vectors
    TV = repmat(rand_coeff_vect(1,:),[Ndim 1]).*repmat(u,[1 size(rand_coeff_vect,2)]) + ...
         repmat(rand_coeff_vect(2,:),[Ndim 1]).*repmat(v,[1 size(rand_coeff_vect,2)]);
    
    for k = 1:size(TV,2)
        
        % translation vector
        tv = TV(:,k);
        I(k,:) = (P1 + tv)';
        
    end
    
    % Add the triangle vertices
    I(end-2,:) = P1';
    I(end-1,:) = P2';
    I(end,:)   = P3';
    I = unique(I','rows')';
    
end


end % sample_and_curve_triangle


%% inflate_triangle_sample_from_opposite_vertex subfunction
function P = inflate_triangle_sample_from_opposite_vertex(U, X, Rho)
%
% Authors : nicolas.douillet (at) free.fr, 2017-2024.
%           Gerd Wachsmuth,                     2021.


% We are looking for t such that
%   || t U - X ||^2 = Rho^2
% This is a quadratic equation.

% Discriminant
D = sum(U.*X, 2).^2 - sum(U.^2,2).*(sum(X.^2,2) - Rho^2);

% We take the positive solution
t = (sum(U.*X, 2) + sqrt(D)) ./ sum(U.^2,2);

P = t.*U;


end % inflate_triangle_sample_from_opposite_vertex