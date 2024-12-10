function V_out = smooth_point_set(V_in, k, nb_it)
%% smooth_point_set : function to smooth one given point set (V_in).
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input arguments
%
%          [| | |]
% - V_in = [X Y Z], real matrix double, the point set, size(V_in) = [nb_vertices,3].
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
% - V_out = [X Y Z], real matrix double, the point set, size(V_out) = [nb_vertices,3].
%           [| | |]


%% Body
V_out = V_in;
nb_vtx = size(V_in,1);


for j = 1:nb_it
    
    % Compute normals
    N = estimate_point_set_normals(V_out,k,'norm','oriented');
    
    % Loop on every point of the set
    for i = 1:nb_vtx
        
        % Search for its k nearest neighbors
        cur_vtx = V_in(i,:);
        ids = knnsearch(V_in,cur_vtx,'k',k,'Distance','seuclidean');
                
        % Compute their average
        G = mean(V_in(ids(1,2:end),:),1);
        N(i,:) = mean(N(ids(1,2:end),:),1);
        
        % Project the current point on the plane orthogonal
        % to its normal and going throught the middle point
        [~,H] = point_to_plane_distance(cur_vtx,N(i,:),G);
        V_out(i,:) = H;
        
    end
        
end


end % smooth_point_set