function P = biscuit(L, l, e, nb_samples, isotropic_sampling, random_sampling)
%% biscuit : function to compute a point set in the shape of a biscuit.
%
%%% Author : nicolas.douillet9 (at)gmail.com, 2016-2025. 
%
%
%%% Input arguments
%
% - L (length, double)
% - l (width, double)
% - e (thickness = 2*sphere_radius, double)
% - nb_samples (integer)
% - isotropic_sampling (bool)
% - random_sampling
%
%
%%% Output argument
%
% - P (data matrix)
% 
% Notes on construction method : 
% 4 quarters of a sphere for corners : S1, S2, S3, S4
% 4 half cylinders for edges : Cdv, Cgv; Cdh, Cgh
% 2 plans for sides Ps, Pi
% 
% Centered in O(0,0)
% Lays in XOY plan
% thickness e in Z direction


%% Body
angle_step = 2*pi/nb_samples;
linstepx = 0.1; % step in x direction % angle_step
linstepy = 0.1; % step in y direction % angle_step

% X & Y sampling vectors
ux = -0.5*L:linstepx:0.5*L;
uy = -0.5*l:linstepy:0.5*l;

% Sphere
r = 0.5*e;
theta = (0:angle_step:pi)';
phi =    0:angle_step:2*pi;
sz_theta1 = size(theta,1);
sz_phi2   = size(phi,2);

if isotropic_sampling
    
    X = @(u,v)r*sin(u).*cos(v);
    Y = @(u,v)r*sin(u).*sin(v);
    Z = @(u,v)r*cos(u);
    
    [K,~,v] = sphere_homeo_sfc_isotropic_splg(X,Y,Z,[0 pi sz_theta1],[0 2*pi sz_phi2],random_sampling);
    
    Sx = reshape(K(:,1), [sz_theta1 sz_phi2]);
    Sy = reshape(K(:,2), [sz_theta1 sz_phi2]);
    Sz = reshape(K(:,3), [sz_theta1 sz_phi2]);
    
    % Split into 4 ("sort" phi vector)
    % and anti-clockwise translations
    
    f1 = find(v < phi(1,end)/4);
    f2 = find(v >= 0.5*pi & v < 0.5*phi(1,end));
    f3 = find(v >= pi & v < 0.75*phi(1,end));
    f4 = find(v >= 1.5*pi & v < phi(1,end));
    
    S1x = 0.5*L + Sx(f1);
    S1y = 0.5*l + Sy(f1);
    S1z = Sz(f1);
    
    S2x = -0.5*L + Sx(f2);
    S2y = 0.5*l + Sy(f2);
    S2z = Sz(f2);
    
    S3x = -0.5*L + Sx(f3);
    S3y = -0.5*l + Sy(f3);
    S3z = Sz(f3);
    
    S4x = 0.5*L + Sx(f4);
    S4y = -0.5*l + Sy(f4);
    S4z = Sz(f4);
    
    
    % Cylinders
    % "Horizontal" (Ox) cylinder
    Chx = L*(rand(sz_phi2, size(ux,2))-0.5);
    hz_iso_vect = 2*pi*rand(sz_phi2, size(Chx,2));
    Chy = r*sin(hz_iso_vect);
    Chz = r*cos(hz_iso_vect);
    
    % "Vertical" (Oy) cylinder
    Cvy = l*(rand(sz_phi2, size(uy,2))-0.5);
    vt_iso_vect = 2*pi*rand(sz_phi2, size(Cvy,2));
    Cvx = r*sin(vt_iso_vect);
    Cvz = r*cos(vt_iso_vect);
    
    % Split into 2 ("sort" phi vector)
    % and translations
    e1 = find(vt_iso_vect < 0.5*phi(1,end));
    e2 = find(vt_iso_vect >= 0.5*phi(1,end));
    e3 = find(hz_iso_vect < 0.5*phi(1,end));
    e4 = find(hz_iso_vect >= 0.5*phi(1,end));
    
    Cdvx = 0.5*L+Cvx(e1);
    Cdvy = Cvy(e1);
    Cdvz = Cvz(e1);
    
    Cgvx = -0.5*L + Cvx(e2);
    Cgvy = Cvy(e2);
    Cgvz = Cvz(e2);
    
    Chhx = Chx(e3);
    Chhy = 0.5*l + Chy(e3);
    Chhz = Chz(e3);
    
    Cbhx = Chx(e4);
    Cbhy = -0.5*l + Chy(e4);
    Cbhz = Chz(e4);
    
    
    % Plans
    % Lower plan
    Pbx = L*(rand(floor(L/linstepx),floor(l/linstepy))-0.5);
    Pby = l*(rand(floor(L/linstepx),floor(l/linstepy))-0.5);
    Pbz = repmat(-0.5*e,[size(Pbx,1) size(Pbx,2)]);
    % Or [Pbx Pby] = meshgrid(...);
    
    % Upper plan
    Phx = L*(rand(floor(L/linstepx),floor(l/linstepy))-0.5);
    Phy = l*(rand(floor(L/linstepx),floor(l/linstepy))-0.5);
    Phz = repmat(0.5*e,[size(Phx,1) size(Phx,2)]);
    
    
else % if ~isotropic_sampling
    
    Sx = r*sin(theta)*cos(phi);
    Sy = r*sin(theta)*sin(phi);
    Sz = repmat(r*cos(theta),[1 length(phi)]);
    
    % Split into 4 and translations
    % anti-clockwise
    half_smplg = 0.5*nb_samples;
    quarter_smplg = 0.25*nb_samples;
    
    S1x = 0.5*L + Sx(:,1:floor(quarter_smplg)+1);
    S1y = 0.5*l + Sy(:,1:floor(quarter_smplg)+1);
    S1z = Sz(:,1:floor(quarter_smplg)+1);
    
    S2x = -0.5*L + Sx(:,floor(quarter_smplg)+1:floor(half_smplg)+1);
    S2y = 0.5*l + Sy(:,floor(quarter_smplg)+1:floor(half_smplg)+1);
    S2z = Sz(:,floor(quarter_smplg)+1:floor(half_smplg)+1);
    
    S3x = -0.5*L + Sx(:,floor(half_smplg)+1:3*floor(quarter_smplg)+1);
    S3y = -0.5*l + Sy(:,floor(half_smplg)+1:3*floor(quarter_smplg)+1);
    S3z = Sz(:,floor(half_smplg)+1:3*floor(quarter_smplg)+1);
    
    S4x = 0.5*L + Sx(:,3*floor(quarter_smplg)+1:end);
    S4y = -0.5*l + Sy(:,3*floor(quarter_smplg)+1:end);
    S4z = Sz(:,3*floor(quarter_smplg)+1:end);
    
    % Cylinders
    % "Horizontal" (Ox) cylinder
    Chx = repmat(ux,[length(phi) 1]);
    Chy = repmat(r*sin(phi'),[1 size(Chx,2)]);
    Chz = repmat(r*cos(phi'),[1 size(Chx,2)]);
    
    % "Vertical" (Oy) cylinder
    Cvy = repmat(uy,[length(phi) 1]);
    Cvx = repmat(r*sin(phi'),[1 size(Cvy,2)]);
    Cvz = repmat(r*cos(phi'),[1 size(Cvy,2)]);
    
    % Split into 2 and translations
    Cdvx = 0.5*L+Cvx(1:1+floor(half_smplg),:);
    Cdvy = Cvy(1:1+floor(half_smplg),:);
    Cdvz = Cvz(1:1+floor(half_smplg),:);
    
    Cgvx = -0.5*L + Cvx(1+floor(half_smplg):end,:);
    Cgvy = Cvy(1+floor(half_smplg):end,:);
    Cgvz = Cvz(1+floor(half_smplg):end,:);
    
    Chhx = Chx(1:1+floor(half_smplg),:);
    Chhy = 0.5*l + Chy(1:1+floor(half_smplg),:);
    Chhz = Chz(1:1+floor(half_smplg),:);
    
    Cbhx = Chx(1+floor(half_smplg):end,:);
    Cbhy = -0.5*l + Chy(1+floor(half_smplg):end,:);
    Cbhz = Chz(1+floor(half_smplg):end,:);
    
    % Plans
    % Lower plan
    Pbx = repmat(ux,[length(uy) 1]);
    Pby = repmat(uy',[1 length(ux)]);
    Pbz = repmat(-e/2,[size(Pbx,1) size(Pbx,2)]);
    % Or [Pbx Pby] = meshgrid(...);
    
    % Upper plan
    Phx = Pbx;
    Phy = Pby;
    Phz = repmat(e/2,[size(Phx,1) size(Phx,2)]);
    
end


% Put together the 10 point set pieces
X = [reshape(Cgvx,[numel(Cgvx) 1]); reshape(Phx,[numel(Phx) 1]); reshape(Cdvx,[numel(Cdvx) 1]); reshape(Pbx,[numel(Pbx) 1]);   reshape(S3x, [numel(S3x) 1]); ...
     reshape(Chhx,[numel(Chhx) 1]); reshape(S2x,[numel(S2x) 1]); reshape(S1x, [numel(S1x) 1]);  reshape(Cbhx,[numel(Cbhx) 1]); reshape(S4x, [numel(S4x) 1])];

Y = [reshape(Cgvy,[numel(Cgvy) 1]); reshape(Phy,[numel(Phy) 1]); reshape(Cdvy,[numel(Cdvy) 1]); reshape(Pby,[numel(Pby) 1]);   reshape(S3y, [numel(S3y) 1]); ...
     reshape(Chhy,[numel(Chhy) 1]); reshape(S2y,[numel(S2y) 1]); reshape(S1y, [numel(S1y) 1]);  reshape(Cbhy,[numel(Cbhy) 1]); reshape(S4y, [numel(S4y) 1])];

Z = [reshape(Cgvz,[numel(Cgvz) 1]); reshape(Phz,[numel(Phz) 1]); reshape(Cdvz,[numel(Cdvz) 1]); reshape(Pbz,[numel(Pbz) 1]);   reshape(S3z, [numel(S3z) 1]); ...
     reshape(Chhz,[numel(Chhz) 1]); reshape(S2z,[numel(S2z) 1]); reshape(S1z, [numel(S1z) 1]);  reshape(Cbhz,[numel(Cbhz) 1]); reshape(S4z, [numel(S4z) 1])];

P = cat(2,X,Y,Z);
P = unique(P,'rows');


end % biscuit