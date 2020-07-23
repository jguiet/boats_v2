%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Computes export attenuation at depth according to Martin's curve 
%**************************************************************************
function [attenuation] = martin_attenuation(depth,ed,b)

    % Remove land
    depth(find(depth>0)) = NaN; 

    % Ignore regions shallower than ephotic depth
    depth(find(depth>-ed)) = -ed;

    % Attenuation map
    attenuation = (-depth/ed).^(b);
end
