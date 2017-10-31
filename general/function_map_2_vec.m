% function_map_2_vec(array_original)
%-----------------------------------------------------------------------------------------
% convert a 2-dimensional map of sites to a vector of sites
% remove land, high latitude, and open-ocean sites (not in an LME)
%-----------------------------------------------------------------------------------------

function [y1 y2 y3] = function_map_2_vec(array_original)

 load /archive/dcarozza/DATA/mask_notlme_high_nan.mat
 nvec = nansum(nansum(~isnan(mask_notlme_high_nan)));

 nlat = 180;
 nlon = 360;

 array_new = nan(nvec,1);
 array_index = nan(nlat,nlon);

 indlat = nan(nvec,1);
 indlon = nan(nvec,1);

 indarray = 1;

 for indi = 1:nlon
   for indj = 1:nlat
          
     if ( ~isnan(mask_notlme_high_nan(indj,indi)) )
       array_new(indarray) = array_original(indj,indi);
       
       indlat(indarray) = indj;
       indlon(indarray) = indi;
       
%       array_index(indj,indi) = indarray;
       indarray = indarray + 1;
     end
   
   end
 end

 y1 = array_new;
 y2 = indlat;
 y3 = indlon;

end % function

%----------------------------------------------------------------------------------------
% END OF SCRIPT