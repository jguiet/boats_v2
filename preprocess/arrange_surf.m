%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Arrange or compute the surface forcing to be compatible with the inputs 
%**************************************************************************
function surface = arrange_surf(surface,surf_unit,lon,lat)

    % Compute surface if it is not provided
    if isa(surface,'char')
        % ATT JG compute surface of grid cells from lon/lat
        surf_unit='[m^2]';
    end
    
    % Convert surface in [m^2]
    if strcmp(surf_unit,'[m^2]')
        % Nothing to change !
    elseif strcmp(surf_unit,'[km^2]')
        surface=km2_2_m2(surface);
    elseif strcmp(surf_unit,'??')
        disp('Aye, unit not considered yet for:')
        disp(surf_unit)
        disp('Ask Jerome...')
        return
    else 
        disp('Something fishy, check your unit for the surface:')
        disp(surf_unit)
        return
    end
        
end