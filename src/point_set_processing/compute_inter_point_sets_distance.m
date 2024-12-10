function C = compute_inter_point_sets_distance(V1,V2,k)
% compute_inter_point_sets_distance : function to compute the "point to
% point" distance  between two bijective point sets.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
% 
% Input arguments
%
%
%
%
% Output arguments
%
%
%


% - C : score colormap of distance function from V1 to V2


% Normalisation des normales :
%
% en option, à la fin pour la normale au vertex considéré seulement, mais
% pas en intermédiaire.


% Body
nb_vtx = size(V1,1);
C = zeros(nb_vtx,1);

for n = 1:nb_vtx
    
    % knn search
    [idx_nearest,dst_nearest] = knnsearch(V2,V1(n,:),'k',k,'distance','euclidean');
    C(n,1) = mean(dst_nearest,1);
    
end


end % compute_inter_point_sets_distance