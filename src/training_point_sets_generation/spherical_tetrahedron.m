function P = spherical_tetrahedron(Rho, r, nb_samples)
%% spherical_tetrahedron : function to compute and save
% a spherical tetrahedron shaped point set P.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
%%% Input arguments
%
% - Rho : positive real scalar double, the large radius.
%
% - r :   positive real scalar double, the small radius.
%
% - nb_samples : positive integer scalar double, the number of samples.
%
% - random_sampling : logical true/false* | 1/0*, the option to random for random sampling.
%
%
%%% Output argument
%
%       [| | |]
% - P = [X Y Z], real matrix double, the point set, size(P) = [nb_points,3].
%       [| | |]


%% Body
delta = asin(-1/3);

% Tetrahedron summit coordinates
A = Rho*[2*sqrt(2) 0 -1]/3;
% B = Rho*[-sqrt(2) sqrt(6) -1]/3;
% C = Rho*[-sqrt(2) -sqrt(6) -1]/3;
S = Rho*[0 0 1];

% Build 1st SA curved segment (angle  = delta)
Z = @(u,v)Rho*sin(u) + r*sin(u).*cos(v);
X = @(u,v)Rho*cos(u) + r*cos(u).*cos(v);
Y = @(u,v)r*sin(v);

% Sampling parameters
range_u = [delta 0.5*pi floor(((0.5*pi-delta)/pi)*nb_samples)+1];
range_v = [0 2*pi floor(0.25*nb_samples)+1];
SA = torus_homeo_sfc_isotropic_splg(X,Y,Z,range_u,range_v,false);

% Z and Y rotation matrices
Mrz = [-0.5 -0.5*sqrt(3) 0;
       0.5*sqrt(3) -0.5  0;
       0 0 1];

Y_angle = 0.5*pi+delta;

Mry = [cos(Y_angle) 0 -sin(Y_angle);
       0            1  0
       sin(Y_angle) 0  cos(Y_angle)];

% Remove obvious extra junctions points I
A_dir_vect = Mry*[1; 0; 0];
e = [];

for k = 1:size(SA,1)
    
    % ref point
    if k <= floor(0.5*size(SA,1))
        
        ref_pt = A;
        
        % vector to data
        vect_2_data = SA(k,:) - ref_pt;
        
        % Sorting points to keep
        if abs(vect_2_data(1,2)/abs(sum(vect_2_data.*A_dir_vect', 2))) <= tan(pi/3)
            
            e = cat(1,e,k);
            
        end
        
    else % if k > floor(0.5*size(SA,1))
        
        ref_pt = S;
        
        % vector to data
        vect_2_data = SA(k,:) - ref_pt;
        
        % Sorting points to keep
        if abs(vect_2_data(1,2)/vect_2_data(1,1)) <= tan(pi/3)
            
            e = cat(1,e,k);
            
        end
        
    end
    
end

SA = SA(e,:);
SB = (Mrz*SA')';
SC = (Mrz*Mrz*SA')';

% AB & AC curved segments
Roa = [5/6 -0.5/sqrt(3) -sqrt(2)/3;
       0.5/sqrt(3) -0.5 sqrt(2)/sqrt(3);
       -sqrt(2)/3 -sqrt(2)/sqrt(3) -1/3];
   
AB = (Roa*SA')';
BC = (Mrz*AB')';
AC = (Mrz*BC')';

P = cat(1,SA,SB,SC,AB,BC,AC);
P = unique(P,'rows');


end % spherical_tetrahedron