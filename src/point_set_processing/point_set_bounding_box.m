function bbox = point_set_bounding_box(P)
%% point_set_bounding_box : function to compute the bounding box of the point set.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input argument
%        
%       [| | |]
% - P = [X Y Z], real matrix double, the point set, size(P) = [nb_points,3].
%       [| | |]
%
%
%%% Output argument
%
% - bbox : real row vector double, the bounding box, size(bbox) = [1,6],
%          bbox = [xmin xmax ymin ymax zmin zmax].


%% Body
tic;

xmin = min(P(:,1));
xmax = max(P(:,1));
ymin = min(P(:,2));
ymax = max(P(:,2));
zmin = min(P(:,3));
zmax = max(P(:,3));

bbox = [xmin xmax ymin ymax zmin zmax];
fprintf('Point set bounding box computed in %d seconds.\n',toc);


end % point_set_bounding_box