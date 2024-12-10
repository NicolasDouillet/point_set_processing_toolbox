function bbox = point_set_bounding_box(V)
%% point_set_bounding_box : function to compute the bounding box of the point set.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input argument
%        
%       [| | |]
% - V = [X Y Z], real matrix double, the point set, size(V) = [nb_vertices,3].
%       [| | |]
%
%
%%% Output argument
%
% - bbox : real row vector double, the bounding box, size(bbox) = [1,6],
%          bbox = [xmin xmax ymin ymax zmin zmax].


%% Body
tic;

xmin = min(V(:,1));
xmax = max(V(:,1));
ymin = min(V(:,2));
ymax = max(V(:,2));
zmin = min(V(:,3));
zmax = max(V(:,3));

bbox = [xmin xmax ymin ymax zmin zmax];
fprintf('Point set bounding box computed in %d seconds.\n',toc);


end % point_set_bounding_box