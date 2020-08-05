%**************************************************************************************************************
% FUNCTION initialize_domains.m
% Initial conditions
% - Set the initial biomass spectrum dfish
% - Set the initial effort (if needed)
%**************************************************************************************************************
function initial = initialize_domains(boats)
  
  %---------------------------------
  % Aliases for variables category for readability
  MAIN=boats.param.main;
  CONV=boats.param.conversion;
  ENVI=boats.param.environment;
  ECOL=boats.param.ecology;
  FORC=boats.forcing;
  STRU=boats.structure;
  
  switch boats.param.main.sim_init     
  %--------------------------------------------------------------------------------------------------------
  % Initial dfish is set by analytical primary-production regime
  %--------------------------------------------------------------------------------------------------------
  case 'PP'            
  
     disp('initialize biomass assuming primary production growth regime');

     %---------------------------------
     % Set npp and temperature maps for dfish initial state   
     % Use annual averages
     npp          = squeeze(nanmean(FORC.npp,3));                          % mmolC m-2 s-1
     npp_ed       = squeeze(nanmean(FORC.npp_ed,3));                       % mmolC m-3 d-1
     pfb          = squeeze(nanmean(FORC.pfb,3));                          % mmolC m-2 s-1
     temp_phyto   = squeeze(nanmean(FORC.temperature_pel,3));                  % degC
     temp_fish_pel    = squeeze(nanmean(FORC.temperature_pel_K,3));                % degK
     temp_fish_dem    = squeeze(nanmean(FORC.temperature_dem_K,3));                % degK

     %---------------------------------
     % Calculate quantities required for dfish 
     s_over_p   = ( -1.0 + ( 1.0 + 4.0 .* npp_ed ./ (exp(ENVI.kappa_eppley(1).*temp_phyto) .* ...
                    ENVI.Prod_star(1)) ).^0.5) .* 0.5;
     frac_lg_du = s_over_p ./ (1.0 + s_over_p);                            % large fraction of PP as in Dunne et al. (2005)
     mphyto     = (ENVI.mc_phy_l.^frac_lg_du) .* (ENVI.mc_phy_s.^(1.0 - frac_lg_du));
  
     mbentho    = ENVI.mc_benthic;  
     
     if (ECOL.pelagic)&&(ECOL.demersal)
         temp_dep_A_P = exp( (-ENVI.E_activation_A(1)/ENVI.k_Boltzmann) .* (1./temp_fish_pel - 1./ENVI.temp_ref_A));
         A_P          = (ECOL.A00(1)/CONV.spery)*temp_dep_A_P;                        % growth rate of Andersen and Beyer (2013, p. 18)
         mortality0_P = (exp(ECOL.zeta1(1))/3)*A_P;
         temp_dep_A_D = exp( (-ENVI.E_activation_A(2)/ENVI.k_Boltzmann) .* (1./temp_fish_dem - 1./ENVI.temp_ref_A));
         A_D          = (ECOL.A00(2)/CONV.spery)*temp_dep_A_D;                        % growth rate of Andersen and Beyer (2013, p. 18)
         mortality0_D = (exp(ECOL.zeta1(2))/3)*A_D;
     else
         temp_dep_A_P = exp( (-ENVI.E_activation_A(1)/ENVI.k_Boltzmann) .* (1./temp_fish_pel - 1./ENVI.temp_ref_A));
         A_P          = (ECOL.A00(1)/CONV.spery)*temp_dep_A_P;                        % growth rate of Andersen and Beyer (2013, p. 18)
         mortality0_P = (exp(ECOL.zeta1(1))/3)*A_P;
     end

     %---------------------------------
     % Map of trophic scaling ---------
     if (ECOL.pelagic)&&(ECOL.demersal)
         initial.tro_sca = ones(FORC.nlat,FORC.nlon,2);
         % Iron limited trophic scaling for Pelagic fish
         te_no3 = ECOL.te(1) * (1 - (FORC.no3min ./ (ECOL.kfe + FORC.no3min) ) );
         initial.tro_sca(:,:,1) = log10(te_no3) ./ log10(ECOL.ppmr(1)); 
         % Normal trophic scaling for Demersal fish
         initial.tro_sca(:,:,2) = initial.tro_sca(:,:,2)*ECOL.tro_sca(2);
     else
         initial.tro_sca = ones(FORC.nlat,FORC.nlon,1);
         % Iron limited trophic scaling for Pelagic fish
         te_no3 = ECOL.te(1) * (1 - (FORC.no3min ./ (ECOL.kfe + FORC.no3min) ) );
         initial.tro_sca(:,:,1) = log10(te_no3) ./ log10(ECOL.ppmr(1));
     end
     
     %---------------------------------
     % Calculate initial dfish
     if (ECOL.pelagic)&&(ECOL.demersal)
         dfish_P      = (1/ECOL.nfish) * (1 - repmat(initial.tro_sca(:,:,1),[1 1 ECOL.nfish ECOL.nfmass])) .* repmat(npp,[1 1 ECOL.nfish ECOL.nfmass]) ./ ...
         ( repmat(mortality0_P,[1 1 ECOL.nfish ECOL.nfmass]) .* repmat(mphyto.^(initial.tro_sca(:,:,1)),[1 1 ECOL.nfish ECOL.nfmass]) .* ...
         STRU.minf_4d.^(ECOL.h_allo(1) + ECOL.b_allo(1) - 1)) .* STRU.fmass_4d(:,:,1:3,:).^(repmat(initial.tro_sca(:,:,1),[1 1 ECOL.nfish ECOL.nfmass]) + ECOL.h_allo(1) - 1);
         dfish_D      = (1/ECOL.nfish) * (1 - repmat(initial.tro_sca(:,:,2),[1 1 ECOL.nfish ECOL.nfmass])) .* repmat(pfb,[1 1 ECOL.nfish ECOL.nfmass]) ./ ...
         ( repmat(mortality0_D,[1 1 ECOL.nfish ECOL.nfmass]) .* repmat(mbentho.^(initial.tro_sca(1,1,2)),[1 1 ECOL.nfish ECOL.nfmass]) .* ...
         STRU.minf_4d.^(ECOL.h_allo(2) + ECOL.b_allo(2) - 1)) .* STRU.fmass_4d(:,:,4:6,:).^(repmat(initial.tro_sca(:,2),[1 1 ECOL.nfish ECOL.nfmass]) + ECOL.h_allo(2) - 1);         

         %---------------------------------
         % Make non existent cells NaNs
         dfish_P(STRU.mask_notexist_4d) = NaN;
         dfish_D(STRU.mask_notexist_4d) = NaN;
         initial.dfish = cat(3,dfish_P,dfish_D);
     else
         dfish_P      = (1/ECOL.nfish) * (1 - repmat(ECOL.tro_sca(:,:,1),[1 1 ECOL.nfish ECOL.nfmass])) .* repmat(npp,[1 1 ECOL.nfish ECOL.nfmass]) ./ ...
         ( repmat(mortality0_P,[1 1 ECOL.nfish ECOL.nfmass]) .* repmat(mphyto.^(ECOL.tro_sca(:,:,1)),[1 1 ECOL.nfish ECOL.nfmass]) .* ...
         STRU.minf_4d.^(ECOL.h_allo(1) + ECOL.b_allo(1) - 1)) .* STRU.fmass_4d.^(repmat(ECOL.tro_sca(:,:,1),[1 1 ECOL.nfish ECOL.nfmass]) + ECOL.h_allo(1) - 1);      

         %---------------------------------
         % Make non existent cells NaNs
         dfish_P(STRU.mask_notexist_4d) = NaN;
         initial.dfish = dfish_P;
     end
     %---------------------------------
     % Economic harvesting (set effort to zero in each group)
     if strcmp(MAIN.sim_type,'h')
         if (ECOL.pelagic)&&(ECOL.demersal)
            initial.effort = zeros(FORC.nlat,FORC.nlon,2*ECOL.nfish);
         else
            initial.effort = zeros(FORC.nlat,FORC.nlon,ECOL.nfish);
         end
     end

  %--------------------------------------------------------------------------------------------------------
  % Initial dfish is set by a restart file  
  %--------------------------------------------------------------------------------------------------------
  otherwise
     %---------------------------------
     % Restart file
     boats_version = boats.param.main.sim_name;
     outdir     = boats.param.path.outdir;
     path_lname_rest = [outdir '/' 'restart_' boats_version '_nh' boats.param.main.sname_rest];

     %---------------------------------
     % Error if specified restart file lname_rest IS NOT IN in the working directory
     if ~exist([path_lname_rest '.mat'],'file')
       error(['Error: restart file ' path_lname_rest '.mat not found']);
     
     %---------------------------------
     % Initial dfish if specified restart file lname_rest IS IN the working directory
     else
       disp(['loading restart file: ' path_lname_rest '.mat']);
       
       %---------------------------------
       % Load restart file
       tmp = load([path_lname_rest '.mat']);
       restart = tmp.restart;
       
       %---------------------------------
       % Use dfish from restart
       initial.dfish  = restart.dfish;
       
       %---------------------------------
       % Economic harvesting
       if strcmp(MAIN.sim_type,'h')
         % Use effort there is a field named effort in the restart file
         if isfield(restart,'effort')
           initial.effort  = restart.effort;
         else
         % Set effort to zero in each group
	    if (ECOL.pelagic)&&(ECOL.demersal) 
                initial.effort = zeros(FORC.nlat,FORC.nlon,2*ECOL.nfish);
	    else
		initial.effort = zeros(FORC.nlat,FORC.nlon,ECOL.nfish);
	    end
         end              
       end
     end
  end
  
%**************************************************************************************************************
% END OF SCRIPT
