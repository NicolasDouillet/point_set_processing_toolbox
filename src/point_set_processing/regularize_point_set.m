function V_out = regularize_point_set(V_in, k, nb_it)
%% regularize_point_set : function to regularize one given point set (V_in).
%
% Author : nicolas.douillet (at) free.fr, 2024.
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
% - nb_it : positive integer scalar double, the number of iterations to
%           perform.
%
%
%%% Output argument
%
%          [| | |]
% - V_out = [X Y Z], real matrix double, the point set, size(V_out) = [nb_vertices,3].
%          [| | |]


%% Body
nb_vtx = size(V_in,1);
ids = zeros(nb_vtx,k);
dst = zeros(nb_vtx,k);
N = sqrt(sum(V_in.^2,2));
V_out = zeros(size(V_in));


for j = 1:nb_it
    
    % Compute raw normals for each vertex of the set
    for i = 1:nb_vtx
        
        % I Look for k nearest neighbors
        [ids(i,:),dst(i,:)] = knnsearch(V_in,V_in(i,:),'k',k,'Distance','seuclidean');
                
        V_ngb = V_in(ids(i,2:end),:);                
        V_out(i,:) = mean(V_ngb,1);
        
    end
    
end

% Keep norm
V_out = N .* V_out ./ sqrt(sum(V_in.^2,2));


end % regularize_point_set