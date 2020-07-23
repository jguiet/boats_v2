%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Compute particle export through the euphotic layer
%**************************************************************************
function [pep] = particle_export(npp,temp,ed)
    % Convert npp    mmolC m-2 s-1-> mmolC m-2 d-1 
    npp_tmp = npp*24*3600;

    % pe (Dunne et al. 2005, eq 1a)
    pe = -0.0101 * temp + 0.0582 * log(npp_tmp/ed) + 0.419;
    
    % Correct for negative pe_values
    pe(find(pe<0)) = 0;
    
    % pep
    pep = pe.*npp;
end
