% function_vec_2_map(vec_original,indlat,indlon)
%-----------------------------------------------------------------------------------------
% convert a vector of sites into a 2-dimensional map of sites
% inverse operation of function_map_2_vec, which was based on
% removing land, high latitude, and open-ocean sites (not in an LME)
%-----------------------------------------------------------------------------------------

function [y1] = function_vec_2_map(vec_original,indlat,indlon)

 load /archive/dcarozza/DATA/mask_notlme_high_nan.mat
 array_new = nan(size(mask_notlme_high_nan));

 for ind = 1:length(indlat)

   array_new(indlat(ind),indlon(ind)) = vec_original(ind);
   
 end

 y1 = array_new;

end % function

%----------------------------------------------------------------------------------------
% END OF SCRIPT