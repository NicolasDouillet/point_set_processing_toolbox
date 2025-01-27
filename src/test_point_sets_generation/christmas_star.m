function V = christmas_star(nb_samples)
%% christmas_star : function to compute
% and save a 3D Christmas star point set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
%%% Input arguments
%
% - nb_samples : positive integer scalar double, sampling > 2. Optional.
%
%
%%% Output arguments
%
%       [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the point set. Size(V) = [nb_points,3].
%       [|  |  | ]


%% Input parsing and default values
if nargin < 1    
    nb_samples = 32;    
end


%% Body

% Geometric parameters
phi_n = 0.5*(1+sqrt(5));
a = 2/sqrt(phi_n*sqrt(5)); % edge length
e = 0.2; % radial sinusoidal variation amplitude

% Z axis 2*pi/5 rotation matrix 
Mrz = [cos(0.4*pi) -sin(0.4*pi) 0;
       sin(0.4*pi)  cos(0.4*pi) 0;
       0            0           1];

centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));
           
% 1st equilateral triangle
V0 = [0 0 1]';
V1 = [sin(centre_angle) 0 cos(centre_angle)]';
V2 = Mrz*V1;

% Bidirectional (u,v) sampling + compute corresponding squared distances vector
P0 = sample_and_shape_triangle(V0,V1,V2,a,nb_samples,e);

% Replicate / rotation -> upper crown
Pu = P0';

for k = 1:4        
    Pu = cat(1,Pu,(Mrz^k*P0)');     
end

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

P1 = sample_and_shape_triangle(V3,V4,V5,a,nb_samples,e);

% Lower crown
Pl = P1';

for k = 1:4      
    Pl = cat(1,Pl,(Mrz^k*P1)');     
end

P = cat(1,Pu,Pl);

% 1st belt triangle
V6 = V1;
V7 = V2;
V8 = Mrz^2*V5;

P2 = sample_and_shape_triangle(V6,V8,V7,a,nb_samples,e);
P2 = P2';

% 2nd belt triangle
V9  = -V6;
V10 = -V7;
V11 = -V8;

P3 = sample_and_shape_triangle(V9,V10,V11,a,nb_samples,e);
P3 = P3';

% Full belt = centre crown
P4 = cat(1,P2,P3);
Pc = P4;

for k = 1:4        
    Pc = cat(1,Pc,(Mrz^k*P4')');      
end

V = cat(1,P,Pc);
V = unique(V,'rows');


end % christmas_star


%% sample_and_shape_triangle subfunction
function V = sample_and_shape_triangle(V1, V2, V3, a, nbstep, e)
%
% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.


% (V1V2, V1V0) basis
u = (V3 - V2);
v = (V1 - V2);
nu = u / norm(u);
nv = v / norm(v);

stepu = a / nbstep;
stepv = a / nbstep;

% Sampling & points creation
V = [];

for m = 0:nbstep
    
    for n = 0:nbstep
        
        if(m+n <= nbstep) % in (V1,V2,V3) triangle conditions ; indices # nb segments                        
       
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;                                          
       
            % Basis (V2V1, V2V3)
            M = V2 + tv;            
            min_dst = min([norm(tv) norm(tv+V2-V1) norm(tv+V2-V3)]);
            
            % Radial vector            
            rv = M + e*cos(2*pi*min_dst/norm(u))*M/norm(M);       
            V = cat(2,V,rv);
            
        end
           
    end

end


end % sample_and_shape_triangle