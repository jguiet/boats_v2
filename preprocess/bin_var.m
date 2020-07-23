%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Binning of a variable at fine resolution on a coarser grid
%**************************************************************************
function [bin_var] = bin_var(var,lon1,lat1,lon2,lat2)

    % Create grid for binning
    vlon2(:,2:size(lon2,2)) = (lon2(:,1:end-1)+lon2(:,2:end))/2;
    vlon2(:,1)              = lon2(:,1)-(lon2(:,2)-lon2(:,1))/2;
    vlon2(:,size(lon2,2)+1) = lon2(:,end)+(lon2(:,end)-lon2(:,end-1))/2;
    vlat2(2:size(lat2,1),:) = (lat2(1:end-1,:)+lat2(2:end,:))/2;
    vlat2(1,:)              = lat2(1,:)-(lat2(2,:)-lat2(1,:))/2;
    vlat2(size(lat2,1)+1,:) = lat2(end,:)+(lat2(end,:)-lat2(end-1,:))/2;

    % Average on coarse grid
    for i = 1:size(lon2,1)
        for j = 1:size(lon2,2)
            % Indexes
            ind_lon = find((lon1>vlon2(i,j))&(lon1<vlon2(i,j+1)));
            ind_lat = find((lat1>vlat2(i,j))&(lat1<vlat2(i+1,j)));
            ind = intersect(ind_lon,ind_lat);
            % Averaging
            bin_var(i,j) = nanmean(var(ind));
	end            
    end
end
