function V_out = remove_outliers(V_in, k, pct2rm)
%% remove_outliers : function to remove rm_pct percentage of the outliers
% present in the given 3D point set V_in.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
% Input arguments
%
% - V_in
%
% - k :      positive integer scalar double, the number of nearest neighbors. k >= size(V_in,1).
%
% - pct2rm : positive integer scalar double in the range |[0; 100]|, the percentage of vertices to remove.
%
%
% Output arguement
%
% - V_out


%% Body
nb_vtx = size(V_in,1);
dst = zeros(nb_vtx,k);

if nb_vtx > 5e4 && license('test','Distrib_Computing_Toolbox')
    
    % delete(gcp('nocreate'));
    new_pool = parpool;
    
    parfor i = 1:nb_vtx
        
        [~,dst(i,:)] = knnsearch(V_in,V_in(i,:),'k',k,'Distance','seuclidean');
        
    end
    
    delete(new_pool);
    
else
    
    for i = 1:nb_vtx
        
        [~,dst(i,:)] = knnsearch(V_in,V_in(i,:),'k',k,'Distance','seuclidean');
        
    end
    
end

knn_mean_dst = mean(dst,2);
[~,vtx_id] = sort(knn_mean_dst,1,'descend');
rm_nb = round(pct2rm*nb_vtx/100);
V_out = V_in(setdiff(1:nb_vtx,vtx_id(1:rm_nb,1)),:);


end % remove_outliers