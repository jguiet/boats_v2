%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Arrange or compute the surface forcing to be compatible with the inputs 
%**************************************************************************
function [npp npp_ed] = arrange_npp(npp,npp_unit,depth,depth_unit,dim,ed,depth_type)

    % Compute npp if it is not provided
    if isa(npp,'char')
        % ATT JG compute npp analytic field if not defined
        npp_unit='[m^2]';
    end
    
    % Compute depth if it is not provided
    if isa(depth,'char')
        % ATT JG compute depth field if not defined
        depth_unit='[m]';
    end    
    
    % Convert npp in [mmolC m^-2 s^-1]
    if strcmp(npp_unit,'[mmolC m^-2 s^-1]')
        % Nothing to change !
    elseif strcmp(npp_unit,'[mmolC m^-2 d^-1]')
        npp=dm1_2_sm1(npp);
    elseif strcmp(npp_unit,'[mgC m^-2 d^-1]')
        npp=mgC_2_mmolC(dm1_2_sm1(npp));
    elseif strcmp(npp_unit,'??')
        disp('Aye, unit not considered yet for:')
        disp(surf_unit)
        disp('Ask Jerome...')
        return
    else 
        disp('Something fishy, check your unit for the npp:')
        disp(npp_unit)
        return
    end
    
    % Convert depth in [m] and prepare
    if strcmp(depth_unit,'[m]')
        % Nothing to change !
    else 
        disp('Something fishy, check your unit for the depth:')
        disp(npp_unit)
        return
    end    
    
    % Prepare depth mask
    if strcmp(depth_type,'map')
        depth(find(depth>75))=ed;
        depth=repmat(depth,[1 1 dim(3)]);
    elseif strcmp(depth_type,'value')
        depth=ones(dim)*ed;
    end
    
    % Distribute npp_ed, the npp over the euphotic zone depth (in  [mmolC m^-2 d^-1])
    npp_ed=npp./depth;
    
    % Convert npp_ed to [mmolC m^-3 d^-1]
    npp_ed=sm1_2_dm1(npp_ed);  
end