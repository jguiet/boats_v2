%**************************************************************************************************************
% FUNCTION set_structure.m
% Define the structure of the domains in various dimensions to save
% computation time
% - Mass class distribution
% - Mass class width
%**************************************************************************************************************
function structure = set_structure(boats,varargin)

 %---------------------------------
 % Aliases for variables category for readability
 ECOL=boats.param.ecology;
 ECON=boats.param.economy;
 FORC=boats.forcing;

 %---------------------------------
 % Mass class structure (bin boundaries in logarithmic space)
 structure.ifmbound     = 1:1:(ECOL.nfmass+1);       	   % size class number bounds: 1,2,3, ...
 structure.fmbound      = ECOL.fmass_0 .* (ECOL.fmass_e/ECOL.fmass_0).^( (structure.ifmbound-1)./(ECOL.nfmass));
 for indm=1:ECOL.nfmass
     structure.fmass(indm) = (structure.fmbound(indm)*structure.fmbound(indm+1))^(0.5);
 end
 % Dummy boundary condition mass class 
 structure.fmass_bc     = structure.fmass(1) * structure.fmass(1)/structure.fmass(2);
 % Width of mass classes (g)
 structure.delfm        = diff(structure.fmbound);				                       
 
 %---------------------------------
 % Save out frequently used arrays and constants repmats of fmass, minf, malpha, delfm, delfm_2end
 % Mass of mass class (g)
 structure.fmass_2d     = repmat(structure.fmass,[ECOL.nfish 1]);
 structure.fmass_4d     = permute(repmat(structure.fmass_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 % Asymptotic mass (g)
 structure.minf_2d      = repmat(ECOL.minf',[1 ECOL.nfmass]);
 structure.minf_4d      = permute(repmat(structure.minf_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 % Maturity mass (g)
 structure.malpha_2d    = repmat(ECOL.malpha',[1 ECOL.nfmass]);
 structure.malpha_4d    = permute(repmat(structure.malpha_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 % Width of mass classes (g)
 structure.delfm_2d     = repmat(structure.delfm,[ECOL.nfish 1]);
 structure.delfm_4d     = permute(repmat(structure.delfm_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 % Width of mass classes from 2nd to final mass class (g)
 structure.delfm_2end_2d = repmat(structure.delfm(2:end),[ECOL.nfish 1]);
 structure.delfm_2end_4d = permute(repmat(structure.delfm_2end_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);

 %---------------------------------
 % Masks
 mask_notexist_2d   = (repmat(structure.fmass,[ECOL.nfish 1]) >  ...
                       repmat(ECOL.minf,[ECOL.nfmass 1])'); 		% NaN where mass is > asymptotic mass
 structure.mask_notexist_4d = permute(repmat(mask_notexist_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 structure.mask_notexist_4d = structure.mask_notexist_4d;                % NaN where mass is > asymptotic mass
 mask_land_2d     = double(FORC.mask(:,:,1));
 mask_land_g_1    = double(repmat(mask_land_2d,[1 1 ECOL.nfish]));        % land mask groups 1 timestep
 mask_land_g_s_1  = double(repmat(mask_land_2d,[1 1 ECOL.nfish ECOL.nfmass])); % land mask integrated mass classes 1 timestep
 structure.mask_land_g_nan  = mask_land_g_1; structure.mask_land_g_nan(mask_land_g_1==1)=nan;
 structure.mask_land_g_s_nan= mask_land_g_s_1; structure.mask_land_g_s_nan(mask_land_g_s_1==1)=nan; 
 
 %-----------------------------------------------------------------------------------------
 % Harvesting selectivity
 structure.selectivity = nan(ECOL.nfish,ECOL.nfmass);
 structure.selectivity(1,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale*ECON.sel_pos_1*ECOL.malpha(1),ECON.sel_slope);
 structure.selectivity(2,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale*ECON.sel_pos_2*ECOL.malpha(2),ECON.sel_slope);
 structure.selectivity(3,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale*ECON.sel_pos_3*ECOL.malpha(3),ECON.sel_slope);
 structure.selectivity_4d   = permute(repmat(structure.selectivity,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 
 %-----------------------------------------------------------------------------------------
 % calculate mass scaling of energy allocation to reproduction  
 rep_scale              = nan(ECOL.nfish,ECOL.nfmass);
 for indf = 1:ECOL.nfish
   rep_scale(indf,:)    = sigmoid_And_mass(structure.fmass,ECOL.rep_pos*ECOL.malpha(indf),ECOL.rep_slope);
 end
 rep_scale_4d           = permute(repmat(rep_scale,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 structure.rep_alloc_frac         = rep_scale_4d .* (1 - ECOL.eff_a) ./ ( (structure.fmass_4d./structure.minf_4d).^(ECOL.b_allo-1) - ECOL.eff_a);
 
 %-----------------------------------------------------------------------------------------
 % Arrays to simplify long calculations
 structure.delfm_4d_over_fmass_4d = structure.delfm_4d ./ structure.fmass_4d;
 structure.minf_4d_p_hplusbm1     = structure.minf_4d.^(ECOL.h_allo + ECOL.b_allo - 1);
 structure.fmass_4d_p_mh          = structure.fmass_4d.^(-ECOL.h_allo);
 structure.minf_4d_p_bm1          = structure.minf_4d.^(ECOL.b_allo-1);
 structure.fmass_4d_p_b           = structure.fmass_4d.^ECOL.b_allo;

%**************************************************************************************************************
% END FUNCTION
