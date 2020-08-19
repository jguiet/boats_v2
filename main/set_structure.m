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
 if (ECOL.pelagic)&&(ECOL.demersal)
     structure.fmass_2d     = repmat(structure.fmass,[2*ECOL.nfish 1]);
     structure.fmass_4d     = permute(repmat(structure.fmass_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 else
     structure.fmass_2d     = repmat(structure.fmass,[ECOL.nfish 1]);
     structure.fmass_4d     = permute(repmat(structure.fmass_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 end
 % Asymptotic mass (g)
 structure.minf_2d      = repmat(ECOL.minf',[1 ECOL.nfmass]);
 structure.minf_4d      = permute(repmat(structure.minf_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 % Maturity mass (g)
 structure.malpha_2d    = repmat(ECOL.malpha',[1 ECOL.nfmass]);
 structure.malpha_4d    = permute(repmat(structure.malpha_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 % Width of mass classes (g)
 if (ECOL.pelagic)&&(ECOL.demersal)
     structure.delfm_2d     = repmat(structure.delfm,[2*ECOL.nfish 1]);
     structure.delfm_4d     = permute(repmat(structure.delfm_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 else
     structure.delfm_2d     = repmat(structure.delfm,[ECOL.nfish 1]);
     structure.delfm_4d     = permute(repmat(structure.delfm_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 end
 % Width of mass classes from 2nd to final mass class (g)
 if (ECOL.pelagic)&&(ECOL.demersal)
     structure.delfm_2end_2d = repmat(structure.delfm(2:end),[2*ECOL.nfish 1]);
     structure.delfm_2end_4d = permute(repmat(structure.delfm_2end_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 else
     structure.delfm_2end_2d = repmat(structure.delfm(2:end),[ECOL.nfish 1]);
     structure.delfm_2end_4d = permute(repmat(structure.delfm_2end_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 end
 %---------------------------------
 % Masks
 mask_notexist_2d   = (repmat(structure.fmass,[ECOL.nfish 1]) >  ...
                       repmat(ECOL.minf,[ECOL.nfmass 1])'); 		% NaN where mass is > asymptotic mass
 structure.mask_notexist_4d = permute(repmat(mask_notexist_2d,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
 mask_land_2d     = double(FORC.mask(:,:,1));
 if (ECOL.pelagic)&&(ECOL.demersal)
     mask_land_g_1    = double(repmat(mask_land_2d,[1 1 2*ECOL.nfish]));        % land mask groups 1 timestep
     mask_land_g_s_1  = double(repmat(mask_land_2d,[1 1 2*ECOL.nfish ECOL.nfmass])); % land mask integrated mass classes 1 timestep
     structure.mask_land_g_nan  = mask_land_g_1; structure.mask_land_g_nan(mask_land_g_1==1)=nan;
     structure.mask_land_g_s_nan= mask_land_g_s_1; structure.mask_land_g_s_nan(mask_land_g_s_1==1)=nan;
     structure.mask_notexist_4d = cat(3,structure.mask_notexist_4d,structure.mask_notexist_4d); % NaN where mass is > asymptotic mass 
 else
     mask_land_g_1    = double(repmat(mask_land_2d,[1 1 ECOL.nfish]));        % land mask groups 1 timestep
     mask_land_g_s_1  = double(repmat(mask_land_2d,[1 1 ECOL.nfish ECOL.nfmass])); % land mask integrated mass classes 1 timestep
     structure.mask_land_g_nan  = mask_land_g_1; structure.mask_land_g_nan(mask_land_g_1==1)=nan;
     structure.mask_land_g_s_nan= mask_land_g_s_1; structure.mask_land_g_s_nan(mask_land_g_s_1==1)=nan; 
 end
 %-----------------------------------------------------------------------------------------
 % Harvesting selectivity
 if (ECOL.pelagic)&&(ECOL.demersal)
     selectivity_P = nan(ECOL.nfish,ECOL.nfmass);
     selectivity_P(1,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(1)*ECON.sel_pos_1*ECOL.malpha(1),ECON.sel_slope(1));
     selectivity_P(2,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(1)*ECON.sel_pos_2*ECOL.malpha(2),ECON.sel_slope(1));
     selectivity_P(3,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(1)*ECON.sel_pos_3*ECOL.malpha(3),ECON.sel_slope(1));
     selectivity_P_4d   = permute(repmat(selectivity_P,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
     selectivity_D = nan(ECOL.nfish,ECOL.nfmass);
     selectivity_D(1,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(2)*ECON.sel_pos_1*ECOL.malpha(1),ECON.sel_slope(2));
     selectivity_D(2,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(2)*ECON.sel_pos_2*ECOL.malpha(2),ECON.sel_slope(2));
     selectivity_D(3,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(2)*ECON.sel_pos_3*ECOL.malpha(3),ECON.sel_slope(2));
     selectivity_D_4d   = permute(repmat(selectivity_D,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);

     structure.selectivity    = cat(1,selectivity_P,selectivity_D);
     structure.selectivity_4d = cat(3,selectivity_P_4d,selectivity_D_4d);
 else
     selectivity_P = nan(ECOL.nfish,ECOL.nfmass);
     selectivity_P(1,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(1)*ECON.sel_pos_1*ECOL.malpha(1),ECON.sel_slope(1));
     selectivity_P(2,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(1)*ECON.sel_pos_2*ECOL.malpha(2),ECON.sel_slope(1));
     selectivity_P(3,:) = sigmoid_And_length(structure.fmass,ECON.sel_pos_scale(1)*ECON.sel_pos_3*ECOL.malpha(3),ECON.sel_slope(1));
     selectivity_P_4d   = permute(repmat(selectivity_P,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
          
     structure.selectivity    = selectivity_P;
     structure.selectivity_4d = selectivity_P_4d;
 end
 
 %-----------------------------------------------------------------------------------------
 % calculate mass scaling of energy allocation to reproduction  
 if (ECOL.pelagic)&&(ECOL.demersal)
     rep_scale_P              = nan(ECOL.nfish,ECOL.nfmass);
     rep_scale_D              = nan(ECOL.nfish,ECOL.nfmass);
     for indf = 1:ECOL.nfish
       rep_scale_P(indf,:)    = sigmoid_And_mass(structure.fmass,ECOL.rep_pos*ECOL.malpha(indf),ECOL.rep_slope);
       rep_scale_D(indf,:)    = sigmoid_And_mass(structure.fmass,ECOL.rep_pos*ECOL.malpha(indf),ECOL.rep_slope);
     end
     rep_scale_P_4d           = permute(repmat(rep_scale_P,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
     rep_scale_D_4d           = permute(repmat(rep_scale_D,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
     rep_alloc_frac_P         = rep_scale_P_4d .* (1 - ECOL.eff_a) ./ ( (structure.fmass_4d(:,:,1:3,:)./structure.minf_4d).^(ECOL.b_allo(1)-1) - ECOL.eff_a);
     rep_alloc_frac_D         = rep_scale_D_4d .* (1 - ECOL.eff_a) ./ ( (structure.fmass_4d(:,:,4:6,:)./structure.minf_4d).^(ECOL.b_allo(2)-1) - ECOL.eff_a);
 
    structure.rep_alloc_frac = cat(3,rep_alloc_frac_P,rep_alloc_frac_D);
 else
     rep_scale_P              = nan(ECOL.nfish,ECOL.nfmass);
     for indf = 1:ECOL.nfish
       rep_scale_P(indf,:)    = sigmoid_And_mass(structure.fmass,ECOL.rep_pos*ECOL.malpha(indf),ECOL.rep_slope);
     end
     rep_scale_P_4d           = permute(repmat(rep_scale_P,[1 1 FORC.nlat FORC.nlon]),[3 4 1 2]);
     rep_alloc_frac_P         = rep_scale_P_4d .* (1 - ECOL.eff_a) ./ ( (structure.fmass_4d./structure.minf_4d).^(ECOL.b_allo(1)-1) - ECOL.eff_a);
 
    structure.rep_alloc_frac = rep_alloc_frac_P;
 end
 
 %-----------------------------------------------------------------------------------------
 % Arrays to simplify long calculations
 structure.delfm_4d_over_fmass_4d = structure.delfm_4d ./ structure.fmass_4d;
 if (ECOL.pelagic)&&(ECOL.demersal)
     minf_4d_p_hplusbm1_P     = structure.minf_4d.^(ECOL.h_allo(1) + ECOL.b_allo(1) - 1);
     fmass_4d_p_mh_P          = structure.fmass_4d(:,:,1:3,:).^(-ECOL.h_allo(1));
     structure.minf_4d_p_bm1_P          = structure.minf_4d.^(ECOL.b_allo(1)-1);
     structure.fmass_4d_p_b_P           = structure.fmass_4d(:,:,1:3,:).^ECOL.b_allo(1);
     minf_4d_p_hplusbm1_D     = structure.minf_4d.^(ECOL.h_allo(2) + ECOL.b_allo(2) - 1);
     fmass_4d_p_mh_D          = structure.fmass_4d(:,:,4:6,:).^(-ECOL.h_allo(2));
     structure.minf_4d_p_bm1_D          = structure.minf_4d.^(ECOL.b_allo(2)-1);
     structure.fmass_4d_p_b_D           = structure.fmass_4d(:,:,4:6,:).^ECOL.b_allo(2);
     
     structure.minf_4d_p_hplusbm1 = cat(3,minf_4d_p_hplusbm1_P,minf_4d_p_hplusbm1_D);
     structure.fmass_4d_p_mh      = cat(3,fmass_4d_p_mh_P,fmass_4d_p_mh_D);
 else
     structure.minf_4d_p_hplusbm1       = structure.minf_4d.^(ECOL.h_allo(1) + ECOL.b_allo(1) - 1);
     structure.fmass_4d_p_mh            = structure.fmass_4d.^(-ECOL.h_allo(1));
     structure.minf_4d_p_bm1_P          = structure.minf_4d.^(ECOL.b_allo(1)-1);
     structure.fmass_4d_p_b_P           = structure.fmass_4d.^ECOL.b_allo(1);
 end

%**************************************************************************************************************
% END FUNCTION
