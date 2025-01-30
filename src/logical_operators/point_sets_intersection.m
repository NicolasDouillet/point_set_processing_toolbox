function [P, i1] = point_sets_intersection(P1, P2, precision)
%% numeric_data_arrays_intersection : function to compute and return the
% intersection between the two point sets P1 and P2.
%
%%% Author : nicolas.douillet9 (at) gmail.com, 2024-2025.


%% Body
if nargin  < 3
   
    precision = 1e3*eps; % default value
    
end

C2 = num2cell(P2,2);

i1 = cellfun(@(r) find(all(abs(r - P1) < precision, 2)), C2,'un',0);
i1 = i1(~cellfun('isempty',i1));
i1 = unique(cell2mat(i1),'stable');

P = P1(i1,:);


end % point_sets_intersection