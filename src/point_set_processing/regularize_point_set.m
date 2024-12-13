function P_out = regularize_point_set(P_in, k, nb_it)
%% regularize_point_set : function to regularize one given point set (P_in).
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
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
% - nb_it : positive integer scalar double, the number of iterations to
%           perform.
%
%
%%% Output argument
%
%          [| | |]
% - P_out = [X Y Z], real matrix double, the point set, size(P_out) = [nb_points,3].
%          [| | |]


%% Body
nb_pt = size(P_in,1);
ids = zeros(nb_pt,k);
dst = zeros(nb_pt,k);
N = sqrt(sum(P_in.^2,2));
P_out = zeros(size(P_in));


for j = 1:nb_it
    
    % Compute raw normals for each point of the set
    for i = 1:nb_pt
        
        % I Look for k nearest neighbors
        [ids(i,:),dst(i,:)] = knnsearch(P_in,P_in(i,:),'k',k,'Distance','seuclidean');
                
        P_ngb = P_in(ids(i,2:end),:);                
        P_out(i,:) = mean(P_ngb,1);
        
    end
    
end

% Keep norm
P_out = N .* P_out ./ sqrt(sum(P_in.^2,2));


end % regularize_point_set