function P = spiky_cell_like_surface(Rho, alpha_l, nb_samples, n, random_sampling)
%% spiky_cell_like_surface : function to compute and save a spiky cell like surface.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2016-2025.
%
%
%%% Input arguments
%
% - Rho : positive real scalar double, the cell radius. Default value is Rho = 1; Optional.
%
% - alpha_l : positive real scalar double, the spikes limit / thickness angle in radian unit. Default value is alpha_l = pi/12; Optional.
%
% - nb_samples : positive integer scalar double, the number of samples. n > 1. Default value is n = 180. Optional.
%
% - n : positive integer scalar double, the number of spike pairs over one round. Default value is n = 4. Optional.
%
% - random_sampling : logical false*/true | 0*/1. Optional.
%
%
%%% Output argument
%
%       [| | |]
% - P = [X Y Z], real matrix double, the point set, size(P) = [nb_points,3].
%       [| | |]


%% Input parsing and default values
if nargin < 5   
    random_sampling = false;    
    if nargin < 4
        n = 4;
        if nargin < 3
            nb_samples = 180;
            if nargin < 2
               alpha_l = pi/12;
               if nargin < 1
                  Rho = 1; 
               end
            end
        end
    end
end


%% Body
% Geometric parameter
r_dst = Rho*alpha_l; % tan(alpha_l)

% Sampling parameters
step = 2*pi / nb_samples;

% rotation angles vectors
theta = (0:step:pi)';
phi = 0:step:2*pi;

if random_sampling
    
    X = @(u,v)Rho*sin(u).*cos(v);
    Y = @(u,v)Rho*sin(u).*sin(v);
    Z = @(u,v)Rho*cos(u);
    
    range_u = [0 pi nb_samples/2+1];
    range_v = [0 2*pi nb_samples];
    
    [P,theta_c,phi_c] = sphere_homeo_sfc_isotropic_splg(X,Y,Z,range_u,range_v,random_sampling);
        
    [P,X_s,Y_s,Z_s,coeff,f] = subfunction(P,Rho,n,r_dst,theta_c,phi_c);            
    
    Rho_s = 0.5*Rho*(1+cos(coeff));
    
    % Option arcos
    X_s = X_s + Rho_s.*sin(theta_c(f)).*cos(phi_c(f));
    Y_s = Y_s + Rho_s.*sin(theta_c(f)).*sin(phi_c(f));
    Z_s = Z_s + Rho_s.*cos(theta_c(f));
    
    P = cat(1,P,cat(2,X_s,Y_s,Z_s));
    
else % if ~random_sampling 
    
    % rotation angles arrays
    theta_c = repmat(theta,[1 length(phi)]);
    phi_c = repmat(phi,[size(theta,1) 1]);
    
    % Spherical coordinates
    X = Rho * sin(theta) * cos(phi);
    Y = Rho * sin(theta) * sin(phi);
    Z = Rho * repmat(cos(theta),[1 length(phi)]);
    
    P = cat(2,X(:),Y(:),Z(:));
        
    [P,X_s,Y_s,Z_s,coeff,f] = subfunction(P,Rho,n,r_dst,theta_c,phi_c);           
    
    % Find and extract top disk / spherical portion
    Delta_z = r_dst*alpha_l;
    Rho_s = 0.5*Rho*(1+cos(coeff));
    
    % Option arcos
    X_s = X_s + Rho_s.*sin(theta_c(f)).*cos(phi_c(f));
    Y_s = Y_s + Rho_s.*sin(theta_c(f)).*sin(phi_c(f));
    Z_s = Z_s + Rho_s.*cos(theta_c(f));
        
    % find top sinusoid
    U = [X_s Y_s Z_s];
    a = find(U(:,3) <= (Rho-Delta_z));
    b = find(abs(U(:,1)) >= r_dst);
    c = find(abs(U(:,2)) >= r_dst);
    g = union(a,b);
    g = union(g,c);
    
    U = U(g,:);
    
    X_s = U(:,1);
    Y_s = U(:,2);
    Z_s = U(:,3);
    
    % Find bottom sinusoid
    V = [X_s Y_s Z_s];
    a = find(V(:,3) >= (Rho-Delta_z));
    b = find(abs(V(:,1)) >= r_dst);
    c = find(abs(V(:,2)) >= r_dst);
    g = union(a,b);
    g = union(g,c);
    
    V = V(g,:);
    
    X_s = V(:,1);
    Y_s = V(:,2);
    Z_s = V(:,3);
    
    % Find an equatorial sinusoid
    W = [X_s Y_s Z_s];
    a = find(W(:,1) >= (Rho-Delta_z));
    b = find(abs(W(:,2)) <= r_dst);
    c = find(abs(W(:,3)) <= r_dst);
    g = intersect(a,b);
    g = intersect(g,c);
    
    S = W(g,:);
    
    % Compute rotations Y (+/-pi/2) and replace
    Rmy = @(angl)[cos(angl) 0 -sin(angl);
                  0         1  0
                  sin(angl) 0  cos(angl)];
    
    Tb = Rmy(0.5*pi) *S'; % top sinusoid
    Bb = Rmy(-0.5*pi)*S'; % bottom sinusoid
    
    P = cat(1,P,W,Tb',Bb');
    
end

P = unique(P,'rows');


end % spiky_cell_like_surface


%% subfunction
function [P, X_s, Y_s, Z_s, coeff, f] = subfunction(P, Rho, n, r_dst,theta_c, phi_c)


theta_c = theta_c(:);
phi_c = phi_c(:);

switch n
    
    case 1
        rotation_angle = pi;
    case 2
        rotation_angle = 0.5*pi;
    case 4
        rotation_angle = 0.25*pi;
    case 8
        rotation_angle = 0.125*pi;
    otherwise
        error('n must be a positive integerin the set {1; 2; 4; 8}');
end

e = [];
f = [];
coeff = [];

for i = 1:size(P,1)
    
    % Closest cylinder basis vertex % floor ? ceil ?
    theta_v = round(100*theta_c(i)/(100*rotation_angle))*rotation_angle;
    phi_v   = round(100*phi_c(i)  /(100*rotation_angle))*rotation_angle;
    
    % Round at the closest rotation_angle multiple vertex
    clst_cyl_bs_vtx = Rho*[sin(theta_v)*cos(phi_v) sin(theta_v)*sin(phi_v) cos(theta_v)];
    
    if norm(P(i,:)-clst_cyl_bs_vtx) >= r_dst
        
        e = cat(1,e,i);
        
    else
        
        f = cat(1,f,i);
        coeff = cat(1,coeff,pi*(norm(P(i,:)-clst_cyl_bs_vtx) / r_dst));
        
    end
    
end

N = P(f,:);
P = P(e,:); % holed version

X_s = N(:,1);
Y_s = N(:,2);
Z_s = N(:,3);


end % subfunction