function N = save_estimate_normals(V, k, mode)
% estimate_normals : function to estimate normals to the points of the given set (V).
%
% Author : nicolas.douillet (at) free.fr, 2024.
%
%
% Input arguments
%
%
% -V
%
%
% - mode : character string in the set : {'raw','norm'*,'RAW','NORM'}, the variable deciding
%          wether to normalize or not the face normals. Case insensitive.
%
%
% Output argument
%
% - N , the vertex normals set.


% Body
nb_vtx = size(V,1);
ids = zeros(nb_vtx,k);
N   = zeros(nb_vtx,3);


% if nb_vtx > 5e4 && license('test','Distrib_Computing_Toolbox')
%
%     % delete(gcp('nocreate'));
%     new_pool = parpool;
%
%     % I Look for k nearest neighbors
%     parfor i = 1:nb_vtx
%
%         ids(i,:) = knnsearch(V,V(i,:),'k',k,'Distance','seuclidean');
%
%     end
%
%     % II Angular sort of k nearest neighbors
%
%
%
%
%     delete(new_pool);
%
% else
    

    % Furthest vertex from isobarycentre
    
    % Orient outward from point cloud
    
    % Neareste meighbor
    
    % Keep the same orientation
    
    % Until the closest vertex to the isobarycentre
    
    % (nearest neighbor path, by removing the previous and the current one each time)


    % I Look for k nearest neighbors
    for i = 1:nb_vtx
        
        ids(i,:) = knnsearch(V,V(i,:),'k',k,'Distance','seuclidean');
                
        % II Angular sort of k nearest neighbors                
        V_ngb = V(ids(i,2:end),:);
        U = V_ngb - V(i,:);
        ref_vect1 = repmat(U(1,:),[k-1,1]);                
        
        bov = cross(U(1,:),U(2,:),2); % No way to control this order
        angl = atan2(sign(dot(cross(ref_vect1,U,2),repmat(bov,[size(U,1),1]),2)).*...
                     sqrt(sum(cross(ref_vect1,U,2).^2,2)),dot(ref_vect1,U,2));
        
        [~,id] = sort(angl);
        S = V_ngb(id,:);
        
        
        % III Link and mesh triangles with current vertex
        T_ngb = cat(2,i*ones(k-1,1),ids(i,2:end)',circshift(ids(i,2:end),-1)');
        
        % Reindex T_ngb according to sorted vertices S using map
        % Triangulations reindexing according to V naturel indices / order
        M = containers.Map(ids(i,:),1:k);
        T_ngb = values(M,num2cell(T_ngb(:)));
        T_ngb = reshape(cell2mat(T_ngb),[k-1,3]);
        
        
        % IV Compute triangle raw normals
        S = cat(1,V(i,:),S);
        Ni = cross(S(T_ngb(:,2),:)-S(T_ngb(:,1),:),S(T_ngb(:,3),:)-S(T_ngb(:,1),:),2);
        
        
        % V Compute vertex normal as the average of these last ones
        N(i,:) = mean(Ni,1);                
        
    end                 
    
    if strcmpi(mode,'norm')
        
        N = N ./ sqrt(sum(N.^2,2));
        
    end
        
    % VI Reorient vertex normals coherently
    G = mean(V,1);
    dst = sqrt(sum(V.^2,2));
    cur_id = find(dst == max(dst),1); % furthest vertex  
    id_vect = cur_id;
    cur_vtx = V(cur_id,:);    
    
    % Normal orientation
    norm_or = sign(dot(N(cur_id,:),cur_vtx - G,2));
    
    if norm_or < 0
        
        N(cur_id,:) = norm_or * N(cur_id,:);        
        
    end            
    
    while  cur_id % ~isequal(sort(id_vect),1:nb_vtx) % n < nb_vtx % ~isequal(numel(idx),nb_vtx) %
                        
        nn_id = knnsearch(V,cur_vtx,'k',3,'Distance','seuclidean');
        
        % First different from the current one         
        nn_id = setdiff(nn_id,id_vect,'stable');
        
        if nn_id
            
            nn_id = nn_id(1);            
            
            if dot(N(nn_id,:),N(cur_id,:),2) < 0 % ~isequal(nxt_norm_or,cur_norm_or) % nxt_norm_or < 0                                
                
                N(nn_id,:) = - N(nn_id,:);
                
            end
            
            cur_id = nn_id;
            cur_vtx = V(cur_id,:);
            V(cur_id,:) = 7*ones(1,3);            
            id_vect = cat(2,id_vect,nn_id);  
            
            % Pb de recherche des voisins pertinents et des points isolés
            % Et des modifs de valeurs des V...
            
            % Trouver un moyen de parcourir tous les points une seule fois, sans les
            % modifier, d'un point à son plus proche voisin
            
            % Il faut donc parcourir tous les points manquants à la fin
            % Et orienter leurs normales en fonction de leurs plus proches
            % voisins non modifiés
             
            
        else
            
            % numel(id_vect) % -> ok
            
            break;
            
        end
        
    end

    
% end













% + bilinear interpolation 2 à deux (pour le lissage)




end % estimate_normals