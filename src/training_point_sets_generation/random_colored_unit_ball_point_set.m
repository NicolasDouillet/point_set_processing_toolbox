function [P, C] = random_colored_unit_ball_point_set(n)
%% random_colored_unit_ball_point_set : function to create one random colored unitary ball point set.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.
%
%
%%% Input argument
%
% - n : positive integer scalar double, the number of points / samples in
%       the point set.
%
%
%%% Output arguments
%
%       [| | |]
% - P = [X Y Z], real matrix double, the point set, size(P) = [nb_points,3].
%       [| | |]
%
%       [| | |]
% - C = [R G B], real matrix double, the color set, size(C) = [nb_points,3].
%       [| | |]


%% Body
P = 2*(rand(n,3)-0.5);
R1 = sqrt(sum(P.^2,2));
R2 = rand(n,1);
P = P .* R2 ./ R1;

cmap = 1-jet.^0.5;

R = sqrt(sum(P.^2,2));
cid = ceil(255*R);
C = cmap(cid,:);


end % random_colored_unit_ball_point_set