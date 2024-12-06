function [dst_array, mean_dst] = distance_map_between_two_bijective_point_sets(V1, V2, option_display)
%
% Author : nicolas.douillet (at) free.fr, 2024.


dst_array = vecnorm((V1-V2)',2)';
mean_dst = mean(dst_array,1);


% plot3 avec colormap
if option_display
    
    plot3()
    
end


end