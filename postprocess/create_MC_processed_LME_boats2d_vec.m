% create_MC_processed_LME_boats2d_vec.m
%-----------------------------------------------------------------------------------------
% Loads output from a Monte Carlo _vec suite and processes it
% Creates output structures with fish, harvest, and effort at the LME scale
% Creates a parameter structure 
%-----------------------------------------------------------------------------------------

clear all

%-----------------------------------------------------------------------------------------
% add paths and matstart
 addpath('/home/jguiet/Matlab_toolbox/')
%  matstart

%-----------------------------------------------------------------------------------------
% load files
% load /archive/dcarozza/DATA/LME_mask.mat
% load /archive/dcarozza/DATA/geo_time.mat
 load /home/jguiet/BOATS_VB1.1/Ecology.mat
% load /Users/jguiet/MODELS/BOATS/src/BOATS_VB1.1_GLOB/Ecology.mat

%-----------------------------------------------------------------------------------------
% directory where output (Monte Carlo suite) is located
 OutDir = '/data/results/jguiet/BOATS/MCRun/data';
% OutDir = '/Users/jguiet/MODELS/BOATS/ensembles/data';
 OutName = '';
 
%-----------------------------------------------------------------------------------------
% directory where processed dataset are saved
 SaveDir = '/data/results/jguiet/BOATS/MCRun';
% SaveDir = '/Users/jguiet/MODELS/BOATS/ensembles';


%-----------------------------------------------------------------------------------------
% parameters, variables, and masks
%  spery     = geo_time.spery;
 surf      = Ecology.surface;
 mask_land = Ecology.mask;

%-----------------------------------------------------------------------------------------
% read shape file LME66.shp
%  S = shaperead('/archive/dcarozza/DATA/LME66/LME66.shp');
 load /home/jguiet/BOATS_VB1.1/LME_mask.mat
% load /Users/jguiet/PROJECTS/Harvest_Peak/Data/LME/LME_mask.mat  % ATT 
% number of LMEs
 nLME = 66;%length(S);

%-----------------------------------------------------------------------------------------
% set up structure of file names to load
 dirData=dir([OutDir '/' OutName '/*.mat']);
 FileNames={dirData.name}';
 nFiles = length(FileNames);
% nFiles = 50;
 disp(['Processing ' num2str(nFiles) ' files']);

%-----------------------------------------------------------------------------------------
% declare array
 nyears         = 400 ;
 ngroups        = 3 ;
 spery          = 31104000 ;

% Set variables to process
 MParams 	    = {'E_activation_A', 'kappa_eppley','Prod_star','tau0','tau1','b_allo','E_activation_m','zeta1','h_allo','A00','egg_surv','sel_slope','sel_pos_scale'};
 nParams = length(MParams);
 MOutNames 	= {'boats.output.annual.fish_gi_t','boats.output.annual.fish_gi_g','boats.output.annual.fish_select_g_out',...
               'boats.output.annual.harvest_gi_t','boats.output.annual.harvest_gi_g',...
               'boats.output.annual.effort_gi_t','boats.output.annual.effort_gi_g'};
 
%-----------------------------------------------------------------------------------------
% Set arrays
%-----------------------------------------------------------------------------------------
% Declare output matrices
 parameters_mc = nan(nParams,nFiles);
 fish_mc       = nan(nyears,nLME,ngroups,nFiles); 
 harvest_mc    = nan(nyears,nLME,ngroups,nFiles); 
 effort_mc     = nan(nyears,nLME,ngroups,nFiles); 

%-----------------------------------------------------------------------------------------
% Find indices corresponding to each LME
 LME_mask_v = squeeze(reshape(LME_mask,[1 180*360]));
 surf_v     = squeeze(reshape(surf,[1 180*360]));
 for iLME=1:nLME
    lme_index{iLME}   = find(LME_mask_v == iLME);
    %disp(['indf: ' num2str(iLME) ' / ' num2str(nLME) '- npoints: ' num2str(length(lme_index{iLME})) ]);
    lme_areas{iLME}   = surf_v(lme_index{iLME});
    lme_areas_a{iLME} = repmat(lme_areas{iLME},[nyears 1 ngroups]);
 end
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Load and process Monte Carlo files
%-----------------------------------------------------------------------------------------
  
 disp(['Loading Output Files']);
  
 for indf=1:nFiles

   %--------------------------------------------------------------------------------------
   % load boats output file
   disp(['indf: ' num2str(indf) ' / ' num2str(nFiles) ]);
   load([OutDir '/' FileNames{indf}]);
    
   %--------------------------------------------------------------------------------------
   % store Monte-Carlo parameter values
   for indp=1:length(MParams)
       if ismember(indp,[1,2,3,7])
           parameters_mc(indp,indf) = boats.param.environment.(MParams{indp});
       elseif ismember(indp,[12,13])
           parameters_mc(indp,indf) = boats.param.economy.(MParams{indp});
       else
           parameters_mc(indp,indf) = boats.param.ecology.(MParams{indp});
       end
   end
    
   %-------------------------------------------------------------------------------------
   % Calculate group harvest and biomass in each of the nLME LMEs for each year    

     for iLME = 1:nLME(1)

      temp1 = squeeze(reshape(boats.output.annual.fish_g_out,[nyears 1 180*360 3]));
      temp2 = temp1(:,lme_index{iLME},:) .* lme_areas_a{iLME};  
      fish_mc(:,iLME,:,indf) = squeeze(nansum(temp2,2));

      temp1 = squeeze(reshape(boats.output.annual.harvest_g_out,[nyears 1 180*360 3]));
      temp2 = temp1(:,lme_index{iLME},:) .* lme_areas_a{iLME};  
      harvest_mc(:,iLME,:,indf) = squeeze(nansum(temp2,2))*spery;
      
      temp1 = squeeze(reshape(boats.output.annual.effort_g_out,[nyears 1 180*360 3]));
      temp2 = temp1(:,lme_index{iLME},:) .* lme_areas_a{iLME};  
      effort_mc(:,iLME,:,indf) = squeeze(nansum(temp2,2));

     end % iLME

 end % indf

%-----------------------------------------------------------------------------------------
% add over groups and LMEs
%-----------------------------------------------------------------------------------------

 fish_total_LME_group    = squeeze(nansum(nansum(fish_mc,2),3));
 harvest_total_LME_group = squeeze(nansum(nansum(harvest_mc,2),3));
 effort_total_LME_group  = squeeze(nansum(nansum(effort_mc,2),3));

%-----------------------------------------------------------------------------------------
% save arrays
%-----------------------------------------------------------------------------------------
 save([SaveDir '/parameters_mc.mat' ],'parameters_mc','-v7.3')
 save([SaveDir '/fish_mc.mat' ],'fish_mc','-v7.3')
 save([SaveDir '/harvest_mc.mat' ],'harvest_mc','-v7.3')
 save([SaveDir '/effort_mc.mat' ],'effort_mc','-v7.3')

 save([SaveDir '/fish_total_LME_group.mat' ],'fish_total_LME_group','-v7.3')
 save([SaveDir '/harvest_total_LME_group.mat' ],'harvest_total_LME_group','-v7.3')
 save([SaveDir '/effort_total_LME_group.mat' ],'effort_total_LME_group','-v7.3')

%-----------------------------------------------------------------------------------------
% END OF SCRIPT
