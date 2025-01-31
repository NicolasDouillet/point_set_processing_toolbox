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

C = 0.5*(P3 + P4);
R = norm(P1 - C);
u = P2;

x_C = C(1,1);
z_C = C(3,1);

alpha = acos(4/3/sqrt(2));

if ~random_sampling
    
    theta = linspace(pi/2-alpha,0,nb_samples+1);
    
    x_SA = R*cos(theta) + x_C;
    z_SA = (abs(R^2 - (x_SA - x_C).^2)).^0.5 + z_C;
    y_SA = zeros(1,length(x_SA));
    
    SA = [x_SA;y_SA;z_SA];
    
    Rmz = [-0.5         -0.5*sqrt(3) 0;
            0.5*sqrt(3) -0.5         0;
            0            0           1];
    
    Rm = [u(1,1)^2-0.5*(1-u(1,1)^2) 1.5*u(1,1)*u(2,1)-u(3,1)*sqrt(3)/2 1.5*u(1,1)*u(3,1)+u(2,1)*sqrt(3)/2;
        1.5*u(1,1)*u(2,1)+u(3,1)*sqrt(3)/2 u(2,1)^2-0.5*(1-u(2,1)^2) 1.5*u(2,1)*u(3,1)-u(1,1)*sqrt(3)/2;
        1.5*u(1,1)*u(3,1)-u(2,1)*sqrt(3)/2 1.5*u(2,1)*u(3,1)+u(1,1)*sqrt(3)/2 u(3,1)^2-0.5*(1-u(3,1)^2)];
    
    SB = Rmz * SA;   
    BA = Rm * Rm * SA;
    
    % Rectangular Coons patch
    nu = length(theta);
    nv = nu;
    
    AB = fliplr(BA);
    
    SAB = zeros(nu,nv,3);
    SAB(:,1,:)   = cat(3,cat(3,SA(1,:)',SA(2,:)'),SA(3,:)');
    SAB(end,:,:) = cat(3,cat(3,AB(1,:)',AB(2,:)'),AB(3,:)');
    SAB(:,end,:) = cat(3,cat(3,SB(1,:)',SB(2,:)'),SB(3,:)');
    SAB(1,:,:)   = repmat(cat(3,cat(3,0,0),1),[1 length(theta) 1]);
    
    
    for n = 1:nu
        
        for p = 1:nv
            
            u = n/(nu-1);
            v = p/(nv-1);
            
            SAB(n,p,1) = (1-u)*SAB(1,p,1)+...
                         (1-v)*SAB(n,1,1)+...
                             u*SAB(end,p,1) +...
                             v*SAB(n,end,1) +...
                   (u-1)*(1-v)*SAB(1,1,1)-...
                           u*v*SAB(end,end,1)+...
                       u*(v-1)*SAB(end,1,1)+...
                       v*(u-1)*SAB(1,end,1);
            
            SAB(n,p,2) = (1-u)*SAB(1,p,2)+...
                         (1-v)*SAB(n,1,2)+...
                             u*SAB(end,p,2) +...
                             v*SAB(n,end,2) +...
                   (u-1)*(1-v)*SAB(1,1,2)-...
                           u*v*SAB(end,end,2)+...
                       u*(v-1)*SAB(end,1,2)+...
                       v*(u-1)*SAB(1,end,2);
            
            SAB(n,p,3) = (1-u)*SAB(1,p,3)+...
                         (1-v)*SAB(n,1,3)+...
                             u*SAB(end,p,3) +...
                             v*SAB(n,end,3) +...
                   (u-1)*(1-v)*SAB(1,1,3)-...
                           u*v*SAB(end,end,3)+...
                       u*(v-1)*SAB(end,1,3)+...
                       v*(u-1)*SAB(1,end,3);
            
        end
        
    end
    
    SAB_x = SAB(:,:,1);
    SAB_y = SAB(:,:,2);
    SAB_z = SAB(:,:,3);
    
    SAB_vector = [SAB_x(:)';SAB_y(:)';SAB_z(:)'];
    
end

SBC = Rmz * SAB_vector;
SCA = Rmz * SBC;
ABC = Rm  * SCA; 

P = [SAB_vector';SBC';SCA';ABC'];
P = unique(P,'rows');


end % Reuleaux_tetrahedron