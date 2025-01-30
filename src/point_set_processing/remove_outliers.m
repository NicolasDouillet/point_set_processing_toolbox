function [P_out, C_out, N_out] = remove_outliers(P_in, k, pct2rm, C_in, N_in)
%% remove_outliers : function to remove rm_pct percentage of the outliers
% present in the given 3D point set P_in.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.
%
%
%%% Input arguments
%
%          [ |    |    |  ]
% - P_in = [X_in Y_in Z_in], real matrix double, the input point set, size(P_in) = [nb_input_points,3]. Mandatory argument.
%          [ |    |    |  ]
%
% - k :      positive integer scalar double, the number of nearest neighbors. k <= size(P_in,1). In practice 4 <= k <= 50.
%
% - pct2rm : positive integer scalar double in the range |[0; 100]|, the percentage of points to remove.
%
%          [ |    |    |  ]
% - C_in = [R_in G_in B_in], integer matrix double, the input color set, size(C_in) = [nb_input_points,3]. Optional argument.
%          [ |    |    |  ]
%
%          [  |     |     |  ]
% - N_in = [Nx_in Ny_in Nz_in], real matrix double, the input normals set, size(N_in) = [nb_input_points,3]. Optional argument.
%          [  |     |     |  ]
%
%
%%% Output arguements
%
%           [  |     |     |  ]
% - P_out = [X_out Y_out Z_out], real matrix double, the output point set, size(P_out) = [nb_output_points,3],
%           [  |     |     |  ]
%
%           where nb_output_points = nb_input_points - nb_duplicata. Mandatory argument.
%
%           [ |     |     |   ]
% - C_out = [R_out G_out B_out], integer matrix double, the output color set, size(C_out) = [nb_output_points,3]. Optional argument.
%           [ |     |     |   ]
%
%           [  |      |      |   ]
% - N_out = [Nx_out Ny_out Nz_out], real matrix double, the output normals set, size(N_in) = [nb_output_points,3]. Optional argument.
%           [  |      |      |   ]


%% Body
nb_pts = size(P_in,1);
dst = zeros(nb_pts,k);

if nb_pts > 5e4 && license('test','Distrib_Computing_Toolbox')
    
    % delete(gcp('nocreate'));
    new_pool = parpool;
    
    parfor i = 1:nb_pts
        
        [~,dst(i,:)] = knnsearch(P_in,P_in(i,:),'k',k,'Distance','seuclidean');
        
    end
    
    delete(new_pool);
    
else
    
    for i = 1:nb_pts
        
        [~,dst(i,:)] = knnsearch(P_in,P_in(i,:),'k',k,'Distance','seuclidean');
        
    end
    
end

knn_mean_dst = mean(dst,2);
[~,pt_id] = sort(knn_mean_dst,1,'descend');
rm_nb = round(pct2rm*nb_pts/100);
P_out = P_in(setdiff(1:nb_pts,pt_id(1:rm_nb,1)),:);

if nargin > 3
   
    C_out = C_in(setdiff(1:nb_pts,pt_id(1:rm_nb,1)),:);
    N_out = N_in(setdiff(1:nb_pts,pt_id(1:rm_nb,1)),:);

    
end


end % remove_outliers