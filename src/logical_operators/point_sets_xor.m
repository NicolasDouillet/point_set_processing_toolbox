function [P, id1, id2] = point_sets_xor(P1, P2, precision)
%% point_sets_xor : function to compute and return the
% exclusive union between the two point sets P1 and P2.
%
%%% Author : nicolas.douillet (at) free.fr, 2024.


%% Body
if nargin  < 3
   
    precision = 1e3*eps; % default value
    
end

C1 = num2cell(P1,2);
C2 = num2cell(P2,2);

i1 = cellfun(@(r) find(all(abs(r - P1) < precision, 2)), C2,'un',0);
i1 = i1(~cellfun('isempty',i1));
i1 = unique(cell2mat(i1),'stable');
id1 = setdiff(1:size(P1,1),i1);

i2 = cellfun(@(r) find(all(abs(r - P2) < precision, 2)), C1,'un',0);
i2 = i2(~cellfun('isempty',i2));
i2 = unique(cell2mat(i2),'stable');
id2 = setdiff(1:size(P2,1),i2);

P = union(P1(id1,:),P2(id2,:),'rows');


end % point_sets_xor