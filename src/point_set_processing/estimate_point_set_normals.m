function N = estimate_point_set_normals(V, k, mode) % orientation as an option
%% estimate_point_set_normals : function to estimate normals to the points of the given set (V).
%
%%% Author : nicolas.douillet (at) free.fr, 2024.
%
%
%%% Input arguments
%
%       [| | |]
% - V = [X Y Z], real matrix double, the point set, size(V) = [nb_vertices,3].
%       [| | |]
%
% - k : positive integer scalar double, the number of neighbor for the k
%       nearest neighbor search. In practive chose k >= 4.
%
% - mode : character string in the set : {'raw','norm'*,'RAW','NORM'}, the variable deciding
%          wether to normalize or not the face normals. Case insensitive.
%
%
%%% Output argument
%
%       [ |  |  |]
% - N : [Nx Ny Nz], real matrix double, the vertex normal vectors, size(N) = [nb_vertices,3].
%       [ |  |  |]


if nargin  < 3
   
        mode = 'norm';
    
end


%% Body
[x_min, idx1] = min(V(:,1));
[x_max, idx2] = max(V(:,1));
[y_min, idy1] = min(V(:,2));
[y_max, idy2] = max(V(:,2));
[z_min, idz1] = min(V(:,3));
[z_max, idz2] = max(V(:,3));

bbc = 0.5*[(x_min+x_max) (y_min+y_max) (z_min+z_max)]; % bounding box centre
bb_idx = [idx1 idx2 idy1 idy2 idz1 idz2]; % fursthest point belongs to this bounding boxe

V = V - bbc; % mean(V,1);

nb_vtx = size(V,1);
ids = zeros(nb_vtx,k);
N   = zeros(nb_vtx,3);

% Compute raw normals for each vertex of the set
for i = 1:nb_vtx
    
    % I Look for k nearest neighbors
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
check_orientation = false(nb_vtx,1);
id_vect = [];

% Orient outward bounding box xyz limits neighbor normals
for i = 1:6
    
    cur_vtx = V(bb_idx(i),:);
    knn_id = knnsearch(V,cur_vtx,'k',k,'Distance','seuclidean');
    
    n2swith_id = nonzeros((dot(N(knn_id,:),V(knn_id,:),2) < 0) .* knn_id');    
    n2swith_id = setdiff(n2swith_id,id_vect);
    
    if ~isempty(n2swith_id)
    
        % Orient outward
        N(n2swith_id,:) = - N(n2swith_id,:);                
        
    end
    
    id_vect = cat(2,id_vect,knn_id);
    check_orientation(knn_id) = true;
    
end


% 8 quadrant bisectrice direction maxima 
sign_vol = prod(V,2);
S = sign(V);
C1_id = (S(:,1) > 0) & (S(:,2) > 0) & (S(:,3) > 0);
C2_id = (S(:,1) > 0) & (S(:,2) > 0) & (S(:,3) < 0);
C3_id = (S(:,1) > 0) & (S(:,2) < 0) & (S(:,3) > 0);
C4_id = (S(:,1) > 0) & (S(:,2) < 0) & (S(:,3) < 0);
C5_id = (S(:,1) < 0) & (S(:,2) > 0) & (S(:,3) > 0);
C6_id = (S(:,1) < 0) & (S(:,2) > 0) & (S(:,3) < 0);
C7_id = (S(:,1) < 0) & (S(:,2) < 0) & (S(:,3) > 0);
C8_id = (S(:,1) < 0) & (S(:,2) < 0) & (S(:,3) < 0);

set1 = sign_vol.*C1_id;
set2 = sign_vol.*C2_id;
set3 = sign_vol.*C3_id;
set4 = sign_vol.*C4_id;
set5 = sign_vol.*C5_id;
set6 = sign_vol.*C6_id;
set7 = sign_vol.*C7_id;
set8 = sign_vol.*C8_id;

id1 = find(set1 == max(set1),1);
id2 = find(set2 == max(set2),1);
id3 = find(set3 == max(set3),1);
id4 = find(set4 == max(set4),1);
id5 = find(set5 == max(set5),1);
id6 = find(set6 == max(set6),1);
id7 = find(set7 == max(set7),1);
id8 = find(set8 == max(set8),1);

quad_ids = unique(cat(2,id1,id2,id3,id4,id5,id6,id7,id8));
quad_ids = setdiff(quad_ids,id_vect);


% Orient outward 8 quadrant bisectrice direction maxima neighbor normals
for i = 1:numel(quad_ids)
    
    cur_vtx = V(quad_ids(i),:);
    knn_id = knnsearch(V,cur_vtx,'k',k,'Distance','seuclidean');
    
    n2swith_id = nonzeros((dot(N(knn_id,:),V(knn_id,:),2) < 0) .* knn_id');    
    n2swith_id = setdiff(n2swith_id,id_vect);
    
    if ~isempty(n2swith_id)
    
        % Orient outward
        N(n2swith_id,:) = - N(n2swith_id,:);                
        
    end
    
    id_vect = cat(2,id_vect,knn_id);
    check_orientation(knn_id) = true;
    
end


dst = sqrt(sum(V.^2,2));
cur_id = find(dst == max(dst),1); % furthest vertex
cur_vtx = V(cur_id,:);

% Normal orientation
norm_or = sign(dot(N(cur_id,:),cur_vtx,2));

if norm_or < 0
    
    % Orient outward
    N(cur_id,:) = norm_or * N(cur_id,:);
    
end


% From neighborhood to neighborhood
while ~isempty(knn_id) % && numel(id_vect) < nb_vtx % cur_id && ~isempty(cur_vtx) &&
    
    knn_id = knnsearch(V,cur_vtx,'k',k,'Distance','seuclidean');
    knn_id = setdiff(knn_id,id_vect,'stable');
               
    n2swith_id = nonzeros((dot(N(knn_id,:),repmat(N(cur_id,:),[numel(knn_id),1]),2) < 0) .* knn_id');    
    n2swith_id = setdiff(n2swith_id,id_vect);
    
    if ~isempty(n2swith_id)
    
        % Orient outward
        N(n2swith_id,:) = - N(n2swith_id,:);                
        
    end
    
    id_vect = cat(2,id_vect,knn_id);
    check_orientation(knn_id) = true;
    
    if ~isempty(knn_id)
        
        cur_id = knn_id(end); % fursthest
        cur_vtx = V(cur_id,:);
        
    end
            
end


rm_ids = setdiff(1:nb_vtx,id_vect);
raw_or_ids   = find(check_orientation);
inter_or_ids = (1:nnz(check_orientation))';
M = containers.Map(inter_or_ids,raw_or_ids);

V_ckeck = V(raw_or_ids,:);
j = 1;

% Orient remaining vertex normals
while j < (1 + numel(rm_ids))
    
    cur_id = rm_ids(j);
    cur_vtx = V(cur_id,:);
    
    % Search nearest neighbor among already normal oriented vertex set
    knn_id = knnsearch(V_ckeck,cur_vtx,'k',k,'Distance','seuclidean');
    
    % Retrieve real indices with map
    knn_id = values(M,num2cell(knn_id));
    knn_id_raw = cell2mat(knn_id);        
    nn_id = knn_id_raw(1,1);
    
    % Nearest neighbor which has already been processed / check
    if ~isempty(nn_id)
        
        if dot(N(nn_id,:),N(cur_id,:),2) < 0
            
            N(cur_id,:) = - N(cur_id,:);
            
        end
        
        check_orientation(cur_id) = true;                
        
    end
    
    j = j + 1;
    
end


end % estimate_point_set_normals