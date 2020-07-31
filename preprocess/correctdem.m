%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Correct the temperature experienced by demersal communities 
%**************************************************************************
function [temp_dem] = correctdem(temp_pel,temp_dem,depth,bathy,ed)
   % For the bottom temperature, replace temperatures in shallow waters (>ed)
   % with surface temperatures. 
   ind = find(bathy>-ed)
   temp_dem(ind) = temp_pel(ind);
end
