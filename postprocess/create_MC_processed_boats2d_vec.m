% create_MC_processed_boats2d_vec.m 
%-----------------------------------------------------------------------------------------
% Loads output from a Monte Carlo _vec suite and processes it
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% add paths and matstart
%  matstart
%  addpath('/archive/dcarozza/MATLAB_PATH/')
%  addpath('/archive/dcarozza/boats_gen/')
 addpath('/home/jguiet/Matlab_toolbox/')

%-----------------------------------------------------------------------------------------
% load files
%  load /archive/dcarozza/DATA/LME_mask.mat
%  load /archive/dcarozza/DATA/geo_time.mat
 load /home/jguiet/BOATS_VB1.1/Ecology.mat
 load /data/results/jguiet/BOATS/MCRun/Obs/LME_FG_NH_2006.mat
 LME_FG_NH_2006 = LME_FG_NH;
 load /data/results/jguiet/BOATS/MCRun/Obs/LME_FG_NH_2010.mat
 LME_FG_NH_2010 = LME_FG_NH;

 load //data/results/jguiet/BOATS/MCRun/Obs/RAM_LME.mat

%-----------------------------------------------------------------------------------------
% load processed MC files
 load /data/results/jguiet/BOATS/MCRun/processed/fish_mc.mat
 load /data/results/jguiet/BOATS/MCRun/processed/harvest_mc.mat
 load /data/results/jguiet/BOATS/MCRun/processed/effort_mc.mat

%-----------------------------------------------------------------------------------------
% parameters
 spery                        = 31104000;
 [nyears nLME ngroups nFiles] = size(harvest_mc);
% nFiles                       = 3

%-----------------------------------------------------------------------------------------
% file names
 dirData=dir(['/data/results/jguiet/BOATS/MCRun/data/*.mat']);
 FileNames={dirData.name}';

%-----------------------------------------------------------------------------------------
% Declare arrays for maxh and y100
%-----------------------------------------------------------------------------------------  

%-----------------------------------------------------------------------------------------
% maxh
 fish_S_maxh             = nan(nFiles,1);
 fish_M_maxh     	     = nan(nFiles,1);
 fish_L_maxh             = nan(nFiles,1);

 harvest_S_maxh          = nan(nFiles,1);
 harvest_M_maxh     	 = nan(nFiles,1);
 harvest_L_maxh          = nan(nFiles,1);

 effort_S_maxh           = nan(nFiles,1);
 effort_M_maxh     		 = nan(nFiles,1);
 effort_L_maxh           = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% y100
 fish_S_y100                      = nan(nFiles,1);
 fish_M_y100     				  = nan(nFiles,1);
 fish_L_y100                      = nan(nFiles,1);

 harvest_S_y100                   = nan(nFiles,1);
 harvest_M_y100     			  = nan(nFiles,1);
 harvest_L_y100                   = nan(nFiles,1);

 effort_S_y100                    = nan(nFiles,1);
 effort_M_y100     				  = nan(nFiles,1);
 effort_L_y100                    = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% Declare arrays for three correlations, p-values, residual, bias, and rmse
% For use within the loop over nFiles
%-----------------------------------------------------------------------------------------  

%-----------------------------------------------------------------------------------------
% harA_nh_max_all_2006
 corr_P_harA_nh_max_all_2006          = nan(nFiles,1);
 corr_S_harA_nh_max_all_2006          = nan(nFiles,1);
 pval_P_harA_nh_max_all_2006          = nan(nFiles,1);
 pval_S_harA_nh_max_all_2006          = nan(nFiles,1);
 residual_harA_nh_max_all_2006        = nan(nFiles,nLME);
 bias_harA_nh_max_all_2006            = nan(nFiles,1);
 rmse_harA_nh_max_all_2006            = nan(nFiles,1);
 
%-----------------------------------------------------------------------------------------
% harA_nh_log10_max_all_2006
 corr_P_harA_nh_log10_max_all_2006          = nan(nFiles,1);
 corr_S_harA_nh_log10_max_all_2006          = nan(nFiles,1);
 pval_P_harA_nh_log10_max_all_2006          = nan(nFiles,1);
 pval_S_harA_nh_log10_max_all_2006          = nan(nFiles,1);
 residual_harA_nh_log10_max_all_2006        = nan(nFiles,nLME);
 bias_harA_nh_log10_max_all_2006            = nan(nFiles,1);
 rmse_harA_nh_log10_max_all_2006            = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harSA_nh_max_all_2006
 corr_P_harSA_nh_max_all_2006        = nan(nFiles,1);
 corr_S_harSA_nh_max_all_2006        = nan(nFiles,1);
 pval_P_harSA_nh_max_all_2006        = nan(nFiles,1);
 pval_S_harSA_nh_max_all_2006        = nan(nFiles,1);
 residual_harSA_nh_max_all_2006      = nan(nFiles,nLME);
 bias_harSA_nh_max_all_2006          = nan(nFiles,1);
 rmse_harSA_nh_max_all_2006          = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harMA_nh_max_all_2006
 corr_P_harMA_nh_max_all_2006        = nan(nFiles,1);
 corr_S_harMA_nh_max_all_2006        = nan(nFiles,1);
 pval_P_harMA_nh_max_all_2006        = nan(nFiles,1);
 pval_S_harMA_nh_max_all_2006        = nan(nFiles,1);
 residual_harMA_nh_max_all_2006      = nan(nFiles,nLME);
 bias_harMA_nh_max_all_2006          = nan(nFiles,1);
 rmse_harMA_nh_max_all_2006          = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harLA_nh_max_all_2006
 corr_P_harLA_nh_max_all_2006        = nan(nFiles,1);
 corr_S_harLA_nh_max_all_2006        = nan(nFiles,1);
 pval_P_harLA_nh_max_all_2006        = nan(nFiles,1);
 pval_S_harLA_nh_max_all_2006        = nan(nFiles,1);
 residual_harLA_nh_max_all_2006      = nan(nFiles,nLME);
 bias_harLA_nh_max_all_2006          = nan(nFiles,1);
 rmse_harLA_nh_max_all_2006          = nan(nFiles,1);
 
%-----------------------------------------------------------------------------------------
% harA_nh_max_all_2010
 corr_P_harA_nh_max_all_2010          = nan(nFiles,1);
 corr_S_harA_nh_max_all_2010          = nan(nFiles,1);
 pval_P_harA_nh_max_all_2010          = nan(nFiles,1);
 pval_S_harA_nh_max_all_2010          = nan(nFiles,1);
 residual_harA_nh_max_all_2010        = nan(nFiles,nLME);
 bias_harA_nh_max_all_2010            = nan(nFiles,1);
 rmse_harA_nh_max_all_2010            = nan(nFiles,1);
 
%-----------------------------------------------------------------------------------------
% harA_nh_log10_max_all_2010
 corr_P_harA_nh_log10_max_all_2010          = nan(nFiles,1);
 corr_S_harA_nh_log10_max_all_2010          = nan(nFiles,1);
 pval_P_harA_nh_log10_max_all_2010          = nan(nFiles,1);
 pval_S_harA_nh_log10_max_all_2010          = nan(nFiles,1);
 residual_harA_nh_log10_max_all_2010        = nan(nFiles,nLME);
 bias_harA_nh_log10_max_all_2010            = nan(nFiles,1);
 rmse_harA_nh_log10_max_all_2010            = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harSA_nh_max_all_2010
 corr_P_harSA_nh_max_all_2010        = nan(nFiles,1);
 corr_S_harSA_nh_max_all_2010        = nan(nFiles,1);
 pval_P_harSA_nh_max_all_2010        = nan(nFiles,1);
 pval_S_harSA_nh_max_all_2010        = nan(nFiles,1);
 residual_harSA_nh_max_all_2010      = nan(nFiles,nLME);
 bias_harSA_nh_max_all_2010          = nan(nFiles,1);
 rmse_harSA_nh_max_all_2010          = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harMA_nh_max_all_2010
 corr_P_harMA_nh_max_all_2010        = nan(nFiles,1);
 corr_S_harMA_nh_max_all_2010        = nan(nFiles,1);
 pval_P_harMA_nh_max_all_2010        = nan(nFiles,1);
 pval_S_harMA_nh_max_all_2010        = nan(nFiles,1);
 residual_harMA_nh_max_all_2010      = nan(nFiles,nLME);
 bias_harMA_nh_max_all_2010          = nan(nFiles,1);
 rmse_harMA_nh_max_all_2010          = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harLA_nh_max_all_2010
 corr_P_harLA_nh_max_all_2010        = nan(nFiles,1);
 corr_S_harLA_nh_max_all_2010        = nan(nFiles,1);
 pval_P_harLA_nh_max_all_2010        = nan(nFiles,1);
 pval_S_harLA_nh_max_all_2010        = nan(nFiles,1);
 residual_harLA_nh_max_all_2010      = nan(nFiles,nLME);
 bias_harLA_nh_max_all_2010          = nan(nFiles,1);
 rmse_harLA_nh_max_all_2010          = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% harvest_2_fish
 harvest_2_fish_min               = nan(nFiles,1);
 harvest_2_fish_mean              = nan(nFiles,1);
 harvest_2_fish_max               = nan(nFiles,1);

%-----------------------------------------------------------------------------------------
% masks

% land mask (2d)
 mask_land      = Ecology.mask;
 mask_high_lat  = LME_FG_NH.mask_high_lat';

%-----------------------------------------------------------------------------------------
% Calculate harvest values from 2006 SAUP data
% Consider SAUP harvest values for the largest years in each LME
% Do this by sorting the total harvest of each LME time series in descending order
%-----------------------------------------------------------------------------------------

 [B I]                          = sort(LME_FG_NH_2006.sum_smlo',2,'descend');
 harvest_SAUP2006_nh_max_all    = nanmean(B(:,1:10),2);
 [B I]                          = sort(LME_FG_NH_2006.harvest_small_adj',2,'descend');
 harvest_SAUP2006_S_nh_max_all  = nanmean(B(:,1:10),2);
 [B I]                          = sort(LME_FG_NH_2006.harvest_medium_adj',2,'descend');
 harvest_SAUP2006_M_nh_max_all  = nanmean(B(:,1:10),2);
 [B I]                          = sort(LME_FG_NH_2006.harvest_large_adj',2,'descend');
 harvest_SAUP2006_L_nh_max_all  = nanmean(B(:,1:10),2);

%-----------------------------------------------------------------------------------------
% Calculate harvest values from 2010 SAUP data
% Consider SAUP harvest values for the largest years in each LME
% Do this by sorting the total harvest of each LME time series in descending order
%-----------------------------------------------------------------------------------------

 [B I]                          = sort(LME_FG_NH_2010.sum_smlo',2,'descend');
 harvest_SAUP2010_nh_max_all    = nanmean(B(:,1:10),2);
 [B I]                          = sort(LME_FG_NH_2010.harvest_small_adj',2,'descend');
 harvest_SAUP2010_S_nh_max_all  = nanmean(B(:,1:10),2);
 [B I]                          = sort(LME_FG_NH_2010.harvest_medium_adj',2,'descend');
 harvest_SAUP2010_M_nh_max_all  = nanmean(B(:,1:10),2);
 [B I]                          = sort(LME_FG_NH_2010.harvest_large_adj',2,'descend');
 harvest_SAUP2010_L_nh_max_all  = nanmean(B(:,1:10),2);

%-----------------------------------------------------------------------------------------
% Calculate harvest, biomass, and harvest to biomass ratios from RAM database
%-----------------------------------------------------------------------------------------

 harvest_RAM                     = RAM_LME.harvest;
 fish_RAM                        = RAM_LME.fish;
 harvest_2_fish_RAM              = harvest_RAM ./ fish_RAM;

%-----------------------------------------------------------------------------------------
% Calculate year of maximum harvest
% Calculate biomass, harvest, and effort at the year of maximum harvest
%-----------------------------------------------------------------------------------------

 fish_total_LME_group    = squeeze(nansum(nansum(fish_mc,2),3));
 harvest_total_LME_group = squeeze(nansum(nansum(harvest_mc,2),3));
 effort_total_LME_group  = squeeze(nansum(nansum(effort_mc,2),3));
 
 [B I]     = sort(harvest_total_LME_group(100:400,:),1,'descend');
 
 maxh      = B(1,:);
 year_maxh = I(1,:)+99;

 for ind = 1:nFiles
 
   fish_S_maxh(ind)       = squeeze(nansum(fish_mc(year_maxh(ind),:,1,ind),2))*1e-15;
   fish_M_maxh(ind)       = squeeze(nansum(fish_mc(year_maxh(ind),:,2,ind),2))*1e-15;
   fish_L_maxh(ind)       = squeeze(nansum(fish_mc(year_maxh(ind),:,3,ind),2))*1e-15;

   harvest_S_maxh(ind)    = squeeze(nansum(harvest_mc(year_maxh(ind),:,1,ind),2))*1e-12;
   harvest_M_maxh(ind)    = squeeze(nansum(harvest_mc(year_maxh(ind),:,2,ind),2))*1e-12;
   harvest_L_maxh(ind)    = squeeze(nansum(harvest_mc(year_maxh(ind),:,3,ind),2))*1e-12;

   effort_S_maxh(ind)     = squeeze(nansum(effort_mc(year_maxh(ind),:,1,ind),2))*1e-9;
   effort_M_maxh(ind)     = squeeze(nansum(effort_mc(year_maxh(ind),:,2,ind),2))*1e-9;
   effort_L_maxh(ind)     = squeeze(nansum(effort_mc(year_maxh(ind),:,3,ind),2))*1e-9;

 end

%-----------------------------------------------------------------------------------------
% Calculate biomass, harvest, and effort at year 100

 for ind = 1:nFiles
 
   fish_S_y100(ind)       = squeeze(nansum(fish_mc(100,:,1,ind),2))*1e-15;
   fish_M_y100(ind)       = squeeze(nansum(fish_mc(100,:,2,ind),2))*1e-15;
   fish_L_y100(ind)       = squeeze(nansum(fish_mc(100,:,3,ind),2))*1e-15;

   harvest_S_y100(ind)    = squeeze(nansum(harvest_mc(100,:,1,ind),2))*1e-12;
   harvest_M_y100(ind)    = squeeze(nansum(harvest_mc(100,:,2,ind),2))*1e-12;
   harvest_L_y100(ind)    = squeeze(nansum(harvest_mc(100,:,3,ind),2))*1e-12;

   effort_S_y100(ind)     = squeeze(nansum(effort_mc(100,:,1,ind),2))*1e-9;
   effort_M_y100(ind)     = squeeze(nansum(effort_mc(100,:,2,ind),2))*1e-9;
   effort_L_y100(ind)     = squeeze(nansum(effort_mc(100,:,3,ind),2))*1e-9;

 end

%-----------------------------------------------------------------------------------------
% Sum over groups
% Select maximum from top 10 years of each LME harvest
 harvest_total_group     = squeeze(nansum(harvest_mc,3));
 [B I]                   = sort(harvest_total_group(100:400,:,:),1,'descend');
 max_harvest_total_group = squeeze(nanmean(B(1:10,:,:),1));

 % sum fish over group
 fish_total_group                = squeeze(nansum(fish_mc,3));

 fish_total_group_prime          = fish_total_group(100:400,:,:);

 % calculate the total group fish over the top ten harvest years from above (use I)
 % calculate mean over these 10 years
 % need to add 99 since the sort occurred over years 100 to 400
 
 fish_total_group_at_max_harvest = nan(10,nLME,nFiles);
 
 for indx = 1:nLME
 
   for indy = 1:nFiles

     fish_total_group_at_max_harvest(:,indx,indy) = fish_total_group_prime(I(1:10,indx,indy),indx,indy);

   end
 end

 fish_total_group_at_max_harvest_mean = squeeze(nanmean(fish_total_group_at_max_harvest,1));

 harvest_total_S         = squeeze(harvest_mc(:,:,1,:));
 [B I]                   = sort(harvest_total_S(100:400,:,:),[1],'descend');
 max_harvest_total_S     = squeeze(nanmean(B(1:10,:,:),1));
 
 harvest_total_M         = squeeze(harvest_mc(:,:,2,:));
 [B I]                   = sort(harvest_total_M(100:400,:,:),[1],'descend');
 max_harvest_total_M     = squeeze(nanmean(B(1:10,:,:),1));

 harvest_total_L         = squeeze(harvest_mc(:,:,3,:));
 [B I]                   = sort(harvest_total_L(100:400,:,:),[1],'descend');
 max_harvest_total_L     = squeeze(nanmean(B(1:10,:,:),1));

%-----------------------------------------------------------------------------------------
% Correlations and pvalues 2006

 for ind = 1:nFiles

 [corr_P_harA_nh_max_all_2006(ind) pval_P_harA_nh_max_all_2006(ind) ...
 corr_S_harA_nh_max_all_2006(ind) pval_S_harA_nh_max_all_2006(ind)] = ...
   nancorrJG(harvest_SAUP2006_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_group(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harSA_nh_max_all_2006(ind) pval_P_harSA_nh_max_all_2006(ind) ...
 corr_S_harSA_nh_max_all_2006(ind) pval_S_harSA_nh_max_all_2006(ind)] = ...
   nancorrJG(harvest_SAUP2006_S_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_S(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harMA_nh_max_all_2006(ind) pval_P_harMA_nh_max_all_2006(ind) ...
 corr_S_harMA_nh_max_all_2006(ind) pval_S_harMA_nh_max_all_2006(ind)] = ...
   nancorrJG(harvest_SAUP2006_M_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_M(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harLA_nh_max_all_2006(ind) pval_P_harLA_nh_max_all_2006(ind) ...
 corr_S_harLA_nh_max_all_2006(ind) pval_S_harLA_nh_max_all_2006(ind)] = ...
   nancorrJG(harvest_SAUP2006_L_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_L(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harA_nh_log10_max_all_2006(ind) pval_P_harA_nh_log10_max_all_2006(ind) ...
 corr_S_harA_nh_log10_max_all_2006(ind) pval_S_harA_nh_log10_max_all_2006(ind)] = ...
   nancorrlog10JG(harvest_SAUP2006_nh_max_all ./ LME_FG_NH.surface_LME,squeeze(max_harvest_total_group(:,ind)) ./ LME_FG_NH.surface_LME);

 end
 
%-----------------------------------------------------------------------------------------
% Correlations and pvalues 2010
 
 for ind = 1:nFiles

 [corr_P_harA_nh_max_all_2010(ind) pval_P_harA_nh_max_all_2010(ind) ...
 corr_S_harA_nh_max_all_2010(ind) pval_S_harA_nh_max_all_2010(ind)] = ...
   nancorrJG(harvest_SAUP2010_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_group(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harSA_nh_max_all_2010(ind) pval_P_harSA_nh_max_all_2010(ind) ...
 corr_S_harSA_nh_max_all_2010(ind) pval_S_harSA_nh_max_all_2010(ind)] = ...
   nancorrJG(harvest_SAUP2010_S_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_S(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harMA_nh_max_all_2010(ind) pval_P_harMA_nh_max_all_2010(ind) ...
 corr_S_harMA_nh_max_all_2010(ind) pval_S_harMA_nh_max_all_2010(ind)] = ...
   nancorrJG(harvest_SAUP2010_M_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_M(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harLA_nh_max_all_2010(ind) pval_P_harLA_nh_max_all_2010(ind) ...
 corr_S_harLA_nh_max_all_2010(ind) pval_S_harLA_nh_max_all_2010(ind)] = ...
   nancorrJG(harvest_SAUP2010_L_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_L(:,ind)) ./ LME_FG_NH.surface_LME);

 [corr_P_harA_nh_log10_max_all_2010(ind) pval_P_harA_nh_log10_max_all_2010(ind) ...
 corr_S_harA_nh_log10_max_all_2010(ind)  pval_S_harA_nh_log10_max_all_2010(ind)] = ...
   nancorrlog10JG(harvest_SAUP2010_nh_max_all ./ LME_FG_NH.surface_LME ,squeeze(max_harvest_total_group(:,ind)) ./ LME_FG_NH.surface_LME);

 end
 
%-----------------------------------------------------------------------------------------
% Harvest to Fish Ratio
 
 harvest_2_fish = max_harvest_total_group ./ fish_total_group_at_max_harvest_mean;

%-------------------------------------------------------------------------------------
% Calculate minimum, maximum, and average H/B for the LMEs in RAM
    
 RAM_LME_mask = logical(repmat(RAM_LME.mask,[1 nFiles]));

 harvest_2_fish(~RAM_LME_mask) = NaN;
    
 for ind = 1:nFiles
    
  harvest_2_fish_min(ind)  = nanmin(harvest_2_fish(:,ind));
  harvest_2_fish_max(ind)  = nanmax(harvest_2_fish(:,ind));
  harvest_2_fish_mean(ind)  = nanmean(harvest_2_fish(:,ind));

 end
 
%%%% ATT JG 
%%%% To Do 
%%%% residual, bias, rmse
%%%% Then entire run script
%%%% Plots
 
%-----------------------------------------------------------------------------------------
% Correlations and pvalues log10

 residual_harA_nh_max_all        = nan(nFiles,nLME);
 bias_harA_nh_max_all            = nan(nFiles,1);
 rmse_harA_nh_max_all            = nan(nFiles,1);

     
 %residual_harSA_nh_max_all(indf,:) = harvest_S_LME_nohigh./LME_FG.surface_LME - harvest_SAUP_S_nohigh_max_all./LME_FG.surface_LME;
 %residual_harSA_nh_max_all(indf,~validdataBoth) = NaN;
 %[bias_harSA_nh_max_all(indf) rmse_harSA_nh_max_all(indf)] = function_residuals(xA,yA);

%-----------------------------------------------------------------------------------------
% Save arrays to MOut
%-----------------------------------------------------------------------------------------

  MOut.nruns     = nFiles;
  MOut.nruns_vec = 1:1:MOut.nruns;

%-----------------------------------------------------------------------------------------
% maxh
  MOut.maxh.fish_S    = fish_S_maxh;
  MOut.maxh.fish_M    = fish_M_maxh;
  MOut.maxh.fish_L    = fish_L_maxh;
  MOut.maxh.fish      = fish_S_maxh + fish_M_maxh + fish_L_maxh;
  MOut.maxh.harvest_S = harvest_S_maxh;
  MOut.maxh.harvest_M = harvest_M_maxh;
  MOut.maxh.harvest_L = harvest_L_maxh;
  MOut.maxh.harvest   = harvest_S_maxh + harvest_M_maxh + harvest_L_maxh;
  MOut.maxh.effort_S  = effort_S_maxh;
  MOut.maxh.effort_M  = effort_M_maxh;
  MOut.maxh.effort_L  = effort_L_maxh;
  MOut.maxh.effort    = effort_S_maxh + effort_M_maxh + effort_L_maxh;
  MOut.maxh.year_maxh = year_maxh;

%-----------------------------------------------------------------------------------------
% y100
  MOut.y100.fish_S    = fish_S_y100;
  MOut.y100.fish_M    = fish_M_y100;
  MOut.y100.fish_L    = fish_L_y100;
  MOut.y100.fish      = fish_S_y100 + fish_M_y100 + fish_L_y100;
  MOut.y100.harvest_S = harvest_S_y100;
  MOut.y100.harvest_M = harvest_M_y100;
  MOut.y100.harvest_L = harvest_L_y100;
  MOut.y100.harvest   = harvest_S_y100 + harvest_M_y100 + harvest_L_y100;
  MOut.y100.effort_S  = effort_S_y100;
  MOut.y100.effort_M  = effort_M_y100;
  MOut.y100.effort_L  = effort_L_y100;
  MOut.y100.effort    = effort_S_y100 + effort_L_y100 + effort_L_y100;

%-----------------------------------------------------------------------------------------
% Save correlations, p-values, residual, bias, and rmse to MOut
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% SAUP 2006

%-----------------------------------------------------------------------------------------
% harA_nh_max_all
  MOut.corr2006.P_harA_nh_max_all             = corr_P_harA_nh_max_all_2006;
  MOut.corr2006.P_R2_harA_nh_max_all          = corr_P_harA_nh_max_all_2006.^2;
  MOut.pval2006.P_harA_nh_max_all             = pval_P_harA_nh_max_all_2006;    
  MOut.corr2006.S_harA_nh_max_all             = corr_S_harA_nh_max_all_2006;
  MOut.pval2006.S_harA_nh_max_all             = pval_S_harA_nh_max_all_2006;
 % MOut.residual2006.harA_nh_max_all           = residual_harA_nh_max_all_2006;
 % MOut.bias2006.harA_nh_max_all               = bias_harA_nh_max_all_2006;
 % MOut.rmse2006.harA_nh_max_all               = rmse_harA_nh_max_all_2006;

%-----------------------------------------------------------------------------------------
% harA_nh_log10_max_all
  MOut.corr2006.P_harA_nh_log10_max_all       = corr_P_harA_nh_log10_max_all_2006;
  MOut.corr2006.P_R2_harA_nh_log10_max_all    = corr_P_harA_nh_log10_max_all_2006.^2;
  MOut.pval2006.P_harA_nh_log10_max_all       = pval_P_harA_nh_log10_max_all_2006;    
 % MOut.residual2006.harA_nh_log10_max_all     = residual_harA_nh_log10_max_all_2006;
 % MOut.bias2006.harA_nh_log10_max_all         = bias_harA_nh_log10_max_all_2006;
 % MOut.rmse2006.harA_nh_log10_max_all         = rmse_harA_nh_log10_max_all_2006;

%-----------------------------------------------------------------------------------------
% harSA_nh_max_all
  MOut.corr2006.P_harSA_nh_max_all             = corr_P_harSA_nh_max_all_2006;
  MOut.corr2006.P_R2_harSA_nh_max_all          = corr_P_harSA_nh_max_all_2006.^2;
  MOut.pval2006.P_harSA_nh_max_all             = pval_P_harSA_nh_max_all_2006;    
  MOut.corr2006.S_harSA_nh_max_all             = corr_S_harSA_nh_max_all_2006;
  MOut.pval2006.S_harSA_nh_max_all             = pval_S_harSA_nh_max_all_2006;
 % MOut.residual2006.harSA_nh_max_all           = residual_harSA_nh_max_all_2006;
 % MOut.bias2006.harSA_nh_max_all               = bias_harSA_nh_max_all_2006;
 % MOut.rmse2006.harSA_nh_max_all               = rmse_harSA_nh_max_all_2006;
  
%-----------------------------------------------------------------------------------------
% harMA_nh_max_all
  MOut.corr2006.P_harMA_nh_max_all             = corr_P_harMA_nh_max_all_2006;
  MOut.corr2006.P_R2_harMA_nh_max_all          = corr_P_harMA_nh_max_all_2006.^2;
  MOut.pval2006.P_harMA_nh_max_all             = pval_P_harMA_nh_max_all_2006;    
  MOut.corr2006.S_harMA_nh_max_all             = corr_S_harMA_nh_max_all_2006;
  MOut.pval2006.S_harMA_nh_max_all             = pval_S_harMA_nh_max_all_2006;
%  MOut.residual2006.harMA_nh_max_all           = residual_harMA_nh_max_all_2006;
%  MOut.bias2006.harMA_nh_max_all               = bias_harMA_nh_max_all_2006;
%  MOut.rmse2006.harMA_nh_max_all               = rmse_harMA_nh_max_all_2006;
  
%-----------------------------------------------------------------------------------------
% harLA_nh_max_all
  MOut.corr2006.P_harLA_nh_max_all             = corr_P_harLA_nh_max_all_2006;
  MOut.corr2006.P_R2_harLA_nh_max_all          = corr_P_harLA_nh_max_all_2006.^2;
  MOut.pval2006.P_harLA_nh_max_all             = pval_P_harLA_nh_max_all_2006;    
  MOut.corr2006.S_harLA_nh_max_all             = corr_S_harLA_nh_max_all_2006;
  MOut.pval2006.S_harLA_nh_max_all             = pval_S_harLA_nh_max_all_2006;
%  MOut.residual2006.harLA_nh_max_all           = residual_harLA_nh_max_all_2006;
%  MOut.bias2006.harLA_nh_max_all               = bias_harLA_nh_max_all_2006;
%  MOut.rmse2006.harLA_nh_max_all               = rmse_harLA_nh_max_all_2006;

%-----------------------------------------------------------------------------------------
% SAUP 2010

%-----------------------------------------------------------------------------------------
% harA_nh_max_all
  MOut.corr2010.P_harA_nh_max_all             = corr_P_harA_nh_max_all_2010;
  MOut.corr2010.P_R2_harA_nh_max_all          = corr_P_harA_nh_max_all_2010.^2;
  MOut.pval2010.P_harA_nh_max_all             = pval_P_harA_nh_max_all_2010;    
  MOut.corr2010.S_harA_nh_max_all             = corr_S_harA_nh_max_all_2010;
  MOut.pval2010.S_harA_nh_max_all             = pval_S_harA_nh_max_all_2010;
%  MOut.residual2010.harA_nh_max_all           = residual_harA_nh_max_all_2010;
%  MOut.bias2010.harA_nh_max_all               = bias_harA_nh_max_all_2010;
%  MOut.rmse2010.harA_nh_max_all               = rmse_harA_nh_max_all_2010;

%-----------------------------------------------------------------------------------------
% harA_nh_log10_max_all
  MOut.corr2010.P_harA_nh_log10_max_all       = corr_P_harA_nh_log10_max_all_2010;
  MOut.corr2010.P_R2_harA_nh_log10_max_all    = corr_P_harA_nh_log10_max_all_2010.^2;
  MOut.pval2010.P_harA_nh_log10_max_all       = pval_P_harA_nh_log10_max_all_2010;
%  MOut.residual2010.harA_nh_log10_max_all     = residual_harA_nh_log10_max_all_2010;
%  MOut.bias2010.harA_nh_log10_max_all         = bias_harA_nh_log10_max_all_2010;
%  MOut.rmse2010.harA_nh_log10_max_all         = rmse_harA_nh_log10_max_all_2010;

%-----------------------------------------------------------------------------------------
% harSA_nh_max_all
  MOut.corr2010.P_harSA_nh_max_all             = corr_P_harSA_nh_max_all_2010;
  MOut.corr2010.P_R2_harSA_nh_max_all          = corr_P_harSA_nh_max_all_2010.^2;
  MOut.pval2010.P_harSA_nh_max_all             = pval_P_harSA_nh_max_all_2010;    
  MOut.corr2010.S_harSA_nh_max_all             = corr_S_harSA_nh_max_all_2010;
  MOut.pval2010.S_harSA_nh_max_all             = pval_S_harSA_nh_max_all_2010;
%  MOut.residual2010.harSA_nh_max_all           = residual_harSA_nh_max_all_2010;
%  MOut.bias2010.harSA_nh_max_all               = bias_harSA_nh_max_all_2010;
%  MOut.rmse2010.harSA_nh_max_all               = rmse_harSA_nh_max_all_2010;
  
%-----------------------------------------------------------------------------------------
% harMA_nh_max_all
  MOut.corr2010.P_harMA_nh_max_all             = corr_P_harMA_nh_max_all_2010;
  MOut.corr2010.P_R2_harMA_nh_max_all          = corr_P_harMA_nh_max_all_2010.^2;
  MOut.pval2010.P_harMA_nh_max_all             = pval_P_harMA_nh_max_all_2010;    
  MOut.corr2010.S_harMA_nh_max_all             = corr_S_harMA_nh_max_all_2010;
  MOut.pval2010.S_harMA_nh_max_all             = pval_S_harMA_nh_max_all_2010;
%  MOut.residual2010.harMA_nh_max_all           = residual_harMA_nh_max_all_2010;
%  MOut.bias2010.harMA_nh_max_all               = bias_harMA_nh_max_all_2010;
%  MOut.rmse2010.harMA_nh_max_all               = rmse_harMA_nh_max_all_2010;
  
%-----------------------------------------------------------------------------------------
% harLA_nh_max_all
  MOut.corr2010.P_harLA_nh_max_all             = corr_P_harLA_nh_max_all_2010;
  MOut.corr2010.P_R2_harLA_nh_max_all          = corr_P_harLA_nh_max_all_2010.^2;
  MOut.pval2010.P_harLA_nh_max_all             = pval_P_harLA_nh_max_all_2010;    
  MOut.corr2010.S_harLA_nh_max_all             = corr_S_harLA_nh_max_all_2010;
  MOut.pval2010.S_harLA_nh_max_all             = pval_S_harLA_nh_max_all_2010;
%  MOut.residual2010.harLA_nh_max_all           = residual_harLA_nh_max_all_2010;
%  MOut.bias2010.harLA_nh_max_all               = bias_harLA_nh_max_all_2010;
%  MOut.rmse2010.harLA_nh_max_all               = rmse_harLA_nh_max_all_2010;

%-----------------------------------------------------------------------------------------
% harvest_2_fish
  MOut.h2f.min                                 = harvest_2_fish_min;
  MOut.h2f.mean                                = harvest_2_fish_mean;
  MOut.h2f.max                                 = harvest_2_fish_max;

% This should go into parameters
  MOut.filenames  = FileNames;

%-----------------------------------------------------------------------------------------
% Analysis
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Define target bounds (harvest units are 10^12 g y-1)
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% define scalings to determine upper and lower ranges of targets (+/- 30%)
 scale_total_low  = 0.5;
 scale_total_high = 1.5;
 
 scale_low        = 0.5;
 scale_high       = 1.5;

 corr		  = 0.45;  
%-----------------------------------------------------------------------------------------
% calculate total and group targets
%-----------------------------------------------------------------------------------------
% take the largest 10 years of total LME harvest
% use these same 10 years for each group
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% 2006

 lme_total_2006        = nansum(LME_FG_NH_2006.sum_smlo,2)*1e-12;
 [B I]                 = sort(lme_total_2006,'descend');
 lme_total_peak_2006   = nanmean(B(1:10));
 
 lme_S_total_2006      = nansum(LME_FG_NH_2006.harvest_small_adj,2)*1e-12;
 lme_S_total_peak_2006 = nanmean(lme_S_total_2006(I(1:10)));

 lme_M_total_2006      = nansum(LME_FG_NH_2006.harvest_medium_adj,2)*1e-12;
 lme_M_total_peak_2006 = nanmean(lme_M_total_2006(I(1:10)));

 lme_L_total_2006      = nansum(LME_FG_NH_2006.harvest_large_adj,2)*1e-12;
 lme_L_total_peak_2006 = nanmean(lme_L_total_2006(I(1:10)));
 
%-----------------------------------------------------------------------------------------
% 2010

 lme_total_2010        = nansum(LME_FG_NH_2010.sum_smlo,2)*1e-12;
 [B I]                 = sort(lme_total_2010,'descend');
 lme_total_peak_2010   = nanmean(B(1:10));
 
 lme_S_total_2010      = nansum(LME_FG_NH_2010.harvest_small_adj,2)*1e-12;
 lme_S_total_peak_2010 = nanmean(lme_S_total_2010(I(1:10)));

 lme_M_total_2010      = nansum(LME_FG_NH_2010.harvest_medium_adj,2)*1e-12;
 lme_M_total_peak_2010 = nanmean(lme_M_total_2010(I(1:10)));

 lme_L_total_2010      = nansum(LME_FG_NH_2010.harvest_large_adj,2)*1e-12;
 lme_L_total_peak_2010 = nanmean(lme_L_total_2010(I(1:10)));

%-----------------------------------------------------------------------------------------
% set target2006
 MOut.target2006.harvest_low    = lme_total_peak_2006*scale_total_low;
 MOut.target2006.harvest_high   = lme_total_peak_2006*scale_total_high;
 MOut.target2006.harvest_S_low  = lme_S_total_peak_2006*scale_low;
 MOut.target2006.harvest_S_high = lme_S_total_peak_2006*scale_high;
 MOut.target2006.harvest_M_low  = lme_M_total_peak_2006*scale_low;
 MOut.target2006.harvest_M_high = lme_M_total_peak_2006*scale_high;
 MOut.target2006.harvest_L_low  = lme_L_total_peak_2006*scale_low;
 MOut.target2006.harvest_L_high = lme_L_total_peak_2006*scale_high;
 MOut.target2006.corr		= 0.67;
%-----------------------------------------------------------------------------------------
% set target2010
 MOut.target2010.harvest_low    = lme_total_peak_2010*scale_low;
 MOut.target2010.harvest_high   = lme_total_peak_2010*scale_high;
 MOut.target2010.harvest_S_low  = lme_S_total_peak_2010*scale_low;
 MOut.target2010.harvest_S_high = lme_S_total_peak_2010*scale_high;
 MOut.target2010.harvest_M_low  = lme_M_total_peak_2010*scale_low;
 MOut.target2010.harvest_M_high = lme_M_total_peak_2010*scale_high;
 MOut.target2010.harvest_L_low  = lme_L_total_peak_2010*scale_low;
 MOut.target2010.harvest_L_high = lme_L_total_peak_2010*scale_high;
 MOut.target2010.corr           = 0.67;
%-----------------------------------------------------------------------------------------
% determine simulations that are within the target bounds
%-----------------------------------------------------------------------------------------
 MOut.target                   = MOut.target2006;
 MOut.corr                     = MOut.corr2006;
 MOut.pval                     = MOut.pval2006;
 

 MOut.sat.harvest              = (MOut.maxh.harvest > MOut.target.harvest_low) & (MOut.maxh.harvest < MOut.target.harvest_high);
 
 MOut.sat.harvest_S            = (MOut.maxh.harvest_S > MOut.target.harvest_S_low) & (MOut.maxh.harvest_S < MOut.target.harvest_S_high);
 MOut.sat.harvest_M            = (MOut.maxh.harvest_M > MOut.target.harvest_M_low) & (MOut.maxh.harvest_M < MOut.target.harvest_M_high);
 MOut.sat.harvest_L            = (MOut.maxh.harvest_L >  MOut.target.harvest_L_low) & (MOut.maxh.harvest_L < MOut.target.harvest_L_high);
 MOut.sat.harvest_all          = MOut.sat.harvest & MOut.sat.harvest_S & MOut.sat.harvest_M & MOut.sat.harvest_L;
  
 MOut.sat.P_A_corr             = (MOut.corr.P_harA_nh_max_all > MOut.target.corr);
 MOut.sat.P_A_corr_harvest     = (MOut.sat.P_A_corr & MOut.sat.harvest);
 MOut.sat.P_A_corr_harvest_all = (MOut.sat.P_A_corr & MOut.sat.harvest_all);
 MOut.sat.S_A_corr             = (MOut.corr.S_harA_nh_max_all > MOut.target.corr);
 MOut.sat.S_A_corr_harvest     = (MOut.sat.S_A_corr & MOut.sat.harvest);
 MOut.sat.S_A_corr_harvest_all = (MOut.sat.S_A_corr & MOut.sat.harvest_all);
 %MOut.sat.K_A_corr             = (MOut.corr.K_harA_nh_max_all > MOut.target.corr);
 %MOut.sat.K_A_corr_harvest     = (MOut.sat.K_A_corr & MOut.sat.harvest);
 %MOut.sat.K_A_corr_harvest_all = (MOut.sat.K_A_corr & MOut.sat.harvest_all);

 MOut.sat.harvest_high         = (MOut.maxh.harvest > MOut.target.harvest_high);
 MOut.sat.harvest_low          = (MOut.maxh.harvest < MOut.target.harvest_low);
 MOut.sat.pval_P_ltp05         = (MOut.pval.P_harA_nh_max_all < 0.05);
 MOut.sat.pval_P_ltp01         = (MOut.pval.P_harA_nh_max_all < 0.01);
 %MOut.sat.pval_K_ltp05         = (MOut.pval.K_harA_nh_max_all < 0.05);
 %MOut.sat.pval_K_ltp01         = (MOut.pval.K_harA_nh_max_all < 0.01);
 MOut.sat.pval_S_ltp05         = (MOut.pval.S_harA_nh_max_all < 0.05);
 MOut.sat.pval_S_ltp01         = (MOut.pval.S_harA_nh_max_all < 0.01);
   
%-----------------------------------------------------------------------------------------
% calculate the number of satisfactory runs and save to MOut
 MOut.satsum.harvest           = sum(MOut.sat.harvest);
 MOut.satsum.harvest_S         = sum(MOut.sat.harvest_S);
 MOut.satsum.harvest_M         = sum(MOut.sat.harvest_M);
 MOut.satsum.harvest_L         = sum(MOut.sat.harvest_L);
 MOut.satsum.harvest_all       = sum(MOut.sat.harvest_all);

 MOut.satsum.P_A_corr             = sum(MOut.sat.P_A_corr);
 MOut.satsum.P_A_corr_harvest     = sum(MOut.sat.P_A_corr_harvest);
 MOut.satsum.P_A_corr_harvest_all = sum(MOut.sat.P_A_corr_harvest_all);
 %MOut.satsum.K_A_corr             = sum(MOut.sat.K_A_corr);
 %MOut.satsum.K_A_corr_harvest     = sum(MOut.sat.K_A_corr_harvest);
 %MOut.satsum.K_A_corr_harvest_all = sum(MOut.sat.K_A_corr_harvest_all);
 MOut.satsum.S_A_corr             = sum(MOut.sat.S_A_corr);
 MOut.satsum.S_A_corr_harvest     = sum(MOut.sat.S_A_corr_harvest);
 MOut.satsum.S_A_corr_harvest_all = sum(MOut.sat.S_A_corr_harvest_all);
  
 MOut.satsum.harvest_high     = sum(MOut.sat.harvest_high);
 MOut.satsum.harvest_low      = sum(MOut.sat.harvest_low);
 MOut.satsum.pval_P_ltp05     = sum(MOut.sat.pval_P_ltp05);
 MOut.satsum.pval_P_ltp01     = sum(MOut.sat.pval_P_ltp01);
 %MOut.satsum.pval_K_ltp05     = sum(MOut.sat.pval_K_ltp05);
 %MOut.satsum.pval_K_ltp01     = sum(MOut.sat.pval_K_ltp01);
 MOut.satsum.pval_S_ltp05     = sum(MOut.sat.pval_S_ltp05);
 MOut.satsum.pval_S_ltp01     = sum(MOut.sat.pval_S_ltp01);

%-----------------------------------------------------------------------------------------
% save harvest_SAUP data to MOut
% MOut.data.harvest_SAUP_nohigh_max_all     = harvest_SAUP_nohigh_max_all;
% MOut.data.harvest_SAUP_S_nohigh_max_all  = harvest_SAUP_S_nohigh_max_all;
% MOut.data.harvest_SAUP_M_nohigh_max_all  = harvest_SAUP_M_nohigh_max_all;
% MOut.data.harvest_SAUP_L_nohigh_max_all  = harvest_SAUP_L_nohigh_max_all;

%-----------------------------------------------------------------------------------------
% find the best corr (Pearson) from among satisfactory total harvest
 MOut.best_corr.P_harA_nh_max5all           = max(MOut.corr.P_harA_nh_max_all(MOut.sat.harvest));
 MOut.best_corr.P_R2_harA_nh_max5all        = max(MOut.corr.P_R2_harA_nh_max_all(MOut.sat.harvest));
 MOut.best_corr.P_harSA                     = max(MOut.corr.P_harSA_nh_max_all(MOut.sat.harvest_S));
 MOut.best_corr.P_R2_harSA                  = max(MOut.corr.P_R2_harSA_nh_max_all(MOut.sat.harvest_S));
 MOut.best_corr.P_harMA                     = max(MOut.corr.P_harMA_nh_max_all(MOut.sat.harvest_M));
 MOut.best_corr.P_R2_harMA                  = max(MOut.corr.P_R2_harMA_nh_max_all(MOut.sat.harvest_M));
 MOut.best_corr.P_harLA                     = max(MOut.corr.P_harLA_nh_max_all(MOut.sat.harvest_L));
 MOut.best_corr.P_R2_harLA                  = max(MOut.corr.P_R2_harLA_nh_max_all(MOut.sat.harvest_L));

% MOut.best_corr.P_har_Sfrac                  = max(MOut.corr.P_Sfrac_nh_max_all(MOut.sat.harvest));
% MOut.best_corr.P_har_Mfrac                  = max(MOut.corr.P_Mfrac_nh_max_all(MOut.sat.harvest));

%-----------------------------------------------------------------------------------------
% find the best corr (Spearman) from among satisfactory harvest
 MOut.best_corr.S_harA_nh_max5all             = max(MOut.corr.S_harA_nh_max_all(MOut.sat.harvest));
 MOut.best_corr.S_harSA                       = max(MOut.corr.S_harSA_nh_max_all(MOut.sat.harvest_S));
 MOut.best_corr.S_harMA                       = max(MOut.corr.S_harMA_nh_max_all(MOut.sat.harvest_M));
 MOut.best_corr.S_harLA                       = max(MOut.corr.S_harLA_nh_max_all(MOut.sat.harvest_L));

% MOut.best_corr.S_har_Sfrac                   = max(MOut.corr.S_Sfrac_nh_max_all(MOut.sat.harvest));
% MOut.best_corr.S_har_Mfrac                   = max(MOut.corr.S_Mfrac_nh_max_all(MOut.sat.harvest));

%-----------------------------------------------------------------------------------------
% find the best corr (Kendall) from among satisfactory harvest
% MOut.best_corr.K_harA_nh_max5all             = max(MOut.corr.K_harA_nh_max_all(MOut.sat.harvest));
% MOut.best_corr.K_harSA                       = max(MOut.corr.K_harSA_nh_max_all(MOut.sat.harvest_S));
% MOut.best_corr.K_harMA                       = max(MOut.corr.K_harMA_nh_max_all(MOut.sat.harvest_M));
% MOut.best_corr.K_harLA                       = max(MOut.corr.K_harLA_nh_max_all(MOut.sat.harvest_L));
%
% MOut.best_corr.K_har_Sfrac                   = max(MOut.corr.K_Sfrac_nh_max_all(MOut.sat.harvest));
% MOut.best_corr.K_har_Mfrac                   = max(MOut.corr.K_Mfrac_nh_max_all(MOut.sat.harvest));

%-----------------------------------------------------------------------------------------
% year
 MOut.year                                    = 1:1:nyears;

%-----------------------------------------------------------------------------------------
% save MOut structure
%-----------------------------------------------------------------------------------------
  save /data/results/jguiet/BOATS/MCRun/MOut_MCV3_0418.mat MOut -v7.3

%-----------------------------------------------------------------------------------------
% END OF SCRIPT
