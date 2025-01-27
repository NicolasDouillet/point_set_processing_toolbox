function [M, u, v] = torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, option_random_sampling)
%% torus_homeo_sfc_isotropic_splg : function to isotropically sample a given torus-homeomorphic surface.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2017-2025.
%
%
%%% Syntax
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z);
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v);
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, option_random_sampling);
%
%
%%% Description
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z)
% generates a tricolon vector / [60*120 3] matrix of X, Y, and Z coordinates
% of points sampling the -torus-homeomorphic- surface defined by
% function handles fctn_x, fctn_y, and fctn_z.
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v)
% generates range_u(1,3)*range_v(1,3) samples located in the area
% [min(u), max(u)] x [min(v) max(v)] = [range_u(1,1), range_u(1,2)] x [range_v(1,1) range_v(1,2)]
%
% torus_homeo_sfc_isotropic_splg(fctn_x, fctn_y, fctn_z, range_u, range_v, option_random_sampling)
% randoms the sampling if option_random_sampling = true/1,
% else -option_random_sampling = false/0- sampling is uniform.
%
%
%%% Input arguments
%
% - fctn_x : function handle in x direction, in spherical coordinates, assumed overloaded for vectors and matrices.
%
% - fctn_y : function handle in y direction, in spherical coordinates, assumed overloaded for vectors and matrices.
%
% - fctn_z : function handle in z direction, in spherical coordinates, assumed overloaded for vectors and matrices.
%
% - range_u : real row vector double, u parameter vector of type : [min(u), max(u), spl_u].
%
% - range_v : real row vector double, v parameter vector of type : [min(v), max(v), spl_v].
%
% - option_random_sampling : logical, *true (1) /false (0).
%
%
%%% Output arguments
%
%       [| | |]
% - M = [X Y Z], real matrix double, the point set. Size(M) = [spl_u*spl_v 3].
%       [| | |]
%
% - u : real matrix double, the sampling matrix / grid in u direction. Size(u) = [spl_u,spl_v].
%
% - v : real matrix double, the sampling matrix / grid in v direction. Size(v) = [spl_u,spl_v].


%% Body
u_min = range_u(1,1);
u_max = range_u(1,2);
spl_u = range_u(1,3);
v_min = range_v(1,1);
v_max = range_v(1,2);
spl_v = range_v(1,3);

[~,u_step] = meshgrid(linspace(0,1,spl_v),linspace(0,1,spl_u));
u_rand = rand(spl_u,spl_v);

u = option_random_sampling * u_rand + (1-option_random_sampling) * u_step;
u = u_min + (u_max-u_min)*u;

v_step = meshgrid(linspace(0,1,spl_v),linspace(0,1,spl_u));
v_rand = rand(spl_u,spl_v);
v = option_random_sampling * v_rand + (1-option_random_sampling) * v_step;
v = v_min + (v_max - v_min)*v;

X = fctn_x(u,v)';
Y = fctn_y(u,v)';
Z = fctn_z(u,v)';
M = cat(2,X(:),Y(:),Z(:));


end % torus_homeo_sfc_isotropic_splg