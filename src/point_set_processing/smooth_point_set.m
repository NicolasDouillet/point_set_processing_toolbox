function P_out = smooth_point_set(P_in, k, nb_it)
%% smooth_point_set : function to smooth one given point set (P_in).
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.
%
%
%%% Input arguments
%
%          [| | |]
% - P_in = [X Y Z], real matrix double, the point set, size(P_in) = [nb_points,3].
%          [| | |]
%
% - k : positive integer scalar double, the number of neighbor for the k
%       nearest neighbor search. In practive chose k >= 4.
%
% - nb_it : positive integer scalar double, the number of smoothing
%           iterations to perform.
%
%
%%% Output argument
%
%
%           [| | |]
% - P_out = [X Y Z], real matrix double, the point set, size(P_out) = [nb_points,3].
%           [| | |]


%% Body
P_out = P_in;
nb_pts = size(P_in,1);


for j = 1:nb_it
    
    % Compute normals
    N = estimate_point_set_normals(P_out,k,'norm','oriented');
    
    % Loop on every point of the set
    for i = 1:nb_pts
        
        % Search for its k nearest neighbors
        cur_pt = P_in(i,:);
        ids = knnsearch(P_in,cur_pt,'k',k,'Distance','seuclidean');
                
        % Compute their average
        G = mean(P_in(ids(1,2:end),:),1);
        N(i,:) = mean(N(ids(1,2:end),:),1);
        
        % Project the current point on the plane orthogonal
        % to its normal and going throught the middle point
        [~,H] = point_to_plane_distance(cur_pt,N(i,:),G);
        P_out(i,:) = H;
        
    end
        
end


end % smooth_point_set