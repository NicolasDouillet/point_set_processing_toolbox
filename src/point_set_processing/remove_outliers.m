function P_out = remove_outliers(P_in, k, pct2rm)
%% remove_outliers : function to remove rm_pct percentage of the outliers
% present in the given 3D point set P_in.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input arguments
%
% - P_in
%
% - k :      positive integer scalar double, the number of nearest neighbors. k >= size(P_in,1).
%
% - pct2rm : positive integer scalar double in the range |[0; 100]|, the percentage of points to remove.
%
%
%%% Output arguement
%
% - P_out


%% Body
nb_pt = size(P_in,1);
dst = zeros(nb_pt,k);

if nb_pt > 5e4 && license('test','Distrib_Computing_Toolbox')
    
    % delete(gcp('nocreate'));
    new_pool = parpool;
    
    parfor i = 1:nb_pt
        
        [~,dst(i,:)] = knnsearch(P_in,P_in(i,:),'k',k,'Distance','seuclidean');
        
    end
    
    delete(new_pool);
    
else
    
    for i = 1:nb_pt
        
        [~,dst(i,:)] = knnsearch(P_in,P_in(i,:),'k',k,'Distance','seuclidean');
        
    end
    
end

knn_mean_dst = mean(dst,2);
[~,pt_id] = sort(knn_mean_dst,1,'descend');
rm_nb = round(pct2rm*nb_pt/100);
P_out = P_in(setdiff(1:nb_pt,pt_id(1:rm_nb,1)),:);


end % remove_outliers