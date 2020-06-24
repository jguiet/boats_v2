% plot_boats2d_MC_process_vec.m
%-----------------------------------------------------------------------------------------
% make plots for an MOut file from a Monte Carlo _vec suite of simulations
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% add paths
% addpath('/archive/dcarozza/MATLAB_PATH/');
 addpath('/home/jguiet/Matlab_toolbox/')
%-----------------------------------------------------------------------------------------
% load processed structure of Monte Carlo simulations
 load /data/results/jguiet/BOATS/MCRun/MOut_MCV3_5.mat
 
% call MOut struc
 struc = MOut;
 nruns = struc.nruns;

% define figs array for plots
 clear figs
 fig_count = 1;

% print options
 idoprint                    = 1;
 idoprint_scatter_FH         = 0;
 idoprint_scatter_param_corr = 0;
 idoprint_histogram          = 0;
 idoprint_qqplot             = 0;
 idoprint_timeseries         = 1;

% label numbers on the LME results
%if do_labels==1
%  label_orig = 1:struc.nruns;
%  labels = num2str(label_orig','%d');
%  text(fish, harvest, labels, 'horizontal','left', 'vertical','bottom','FontSize',8)  
%else
%  figure; scatter(log10(x),log10(y),[100],1:length(x),'o','fill');
%  colorbar
%end

%-----------------------------------------------------------------------------------------
% load variables from structure
   
% calculations_maxh
 fish      = struc.maxh.fish;
 fish_S    = struc.maxh.fish_S;
 fish_M    = struc.maxh.fish_M;
 fish_L    = struc.maxh.fish_L;
 harvest   = struc.maxh.harvest;
 harvest_S = struc.maxh.harvest_S;
 harvest_M = struc.maxh.harvest_M;
 harvest_L = struc.maxh.harvest_L;
 effort    = struc.maxh.effort;
 effort_S  = struc.maxh.effort_S;
 effort_M  = struc.maxh.effort_M;
 effort_L  = struc.maxh.effort_L;

%-----------------------------------------------------------------------------------------
% target bounds
 target_harvest_low           = struc.target.harvest_low;
 target_harvest_high          = struc.target.harvest_high;
 target_harvest_S_low         = struc.target.harvest_S_low;
 target_harvest_S_high        = struc.target.harvest_S_high;
 target_harvest_M_low         = struc.target.harvest_M_low;
 target_harvest_M_high        = struc.target.harvest_M_high;
 target_harvest_L_low         = struc.target.harvest_L_low;
 target_harvest_L_high        = struc.target.harvest_L_high;
 target_harvest_corr          = struc.target.corr;

%-----------------------------------------------------------------------------------------
% satisfactory harvests
 harvest_sat          = struc.sat.harvest;
 harvest_S_sat        = struc.sat.harvest_S;
 harvest_M_sat        = struc.sat.harvest_M;
 harvest_L_sat        = struc.sat.harvest_L;
 harvest_all          = struc.sat.harvest_all;
 P_A_corr             = struc.sat.P_A_corr;
 P_A_corr_harvest     = struc.sat.P_A_corr_harvest;
 P_A_corr_harvest_all = struc.sat.P_A_corr_harvest_all;
 S_A_corr             = struc.sat.S_A_corr;
 S_A_corr_harvest     = struc.sat.S_A_corr_harvest;
 S_A_corr_harvest_all = struc.sat.S_A_corr_harvest_all;
 K_A_corr             = struc.sat.K_A_corr;
 K_A_corr_harvest     = struc.sat.K_A_corr_harvest;
 K_A_corr_harvest_all = struc.sat.K_A_corr_harvest_all;
 harvest_high         = struc.sat.harvest_high;
 harvest_low          = struc.sat.harvest_low;
 pval_P_ltp05         = struc.sat.pval_P_ltp05;
 pval_P_ltp01         = struc.sat.pval_P_ltp01;
 pval_K_ltp05         = struc.sat.pval_K_ltp05;
 pval_K_ltp01         = struc.sat.pval_K_ltp01;
 pval_S_ltp05         = struc.sat.pval_S_ltp05;
 pval_S_ltp01         = struc.sat.pval_S_ltp01;
  
%-----------------------------------------------------------------------------------------
% list variables for plots
%-----------------------------------------------------------------------------------------

 ind_vcv = 1;

%-----------------------------------------------------------------------------------------
% parameters

 var_color_vec{ind_vcv}  = struc.E_activation_A;
 var_string_vec{ind_vcv} = 'Activation energy (\omega_{a,A})'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.kappa_eppley;
 var_string_vec{ind_vcv} = 'Eppley constant';                  ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.Prod_star;
 var_string_vec{ind_vcv} = 'Production star';                  ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.te;
 var_string_vec{ind_vcv} = 'Trophic Efficiency';               ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.ppmr;
 var_string_vec{ind_vcv} = 'PPMR';                             ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = log10(struc.te)./log10(struc.ppmr);
 var_string_vec{ind_vcv} = 'Trophic scaling';                  ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.b_allo;
 var_string_vec{ind_vcv} = 'Allometric growth constant';       ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.E_activation_m;
 var_string_vec{ind_vcv} = 'Activation energy (mortality)';    ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.zeta1;
 var_string_vec{ind_vcv} = 'Mortality zeta1';                  ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.h_allo;
 var_string_vec{ind_vcv} = 'Mortality scaling (h_allo)';       ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.A00;
 var_string_vec{ind_vcv} = 'A00';                              ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.egg_surv;
 var_string_vec{ind_vcv} = 'Egg survival';                     ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = log10(struc.egg_surv);
 var_string_vec{ind_vcv} = 'log10 Egg survival';               ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.sel_slope;
 var_string_vec{ind_vcv} = 'Selectivity slope';                ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.sel_pos_scale;
 var_string_vec{ind_vcv} = 'Selectivity position scaling';     ind_vcv = ind_vcv + 1;

%-----------------------------------------------------------------------------------------
% correlations

 var_color_vec{ind_vcv}  = struc.corr.P_R2_harA_nh_max5_all;
 var_string_vec{ind_vcv} = 'PEARSON R2 HAR PER AREA NH MAX5 ALL';  ind_vcv = ind_vcv + 1;
% var_color_vec{ind_vcv}  = struc.corr.K_harA_nh_max5_all;
% var_string_vec{ind_vcv} = 'KENDAL TAU HAR PER AREA NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_harA_nh_max5_all;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO HAR PER AREA NH MAX5 ALL'; ind_vcv = ind_vcv + 1;

%var_color_vec{ind_vcv} = struc.corr.P_harvest_2_fish;
%var_string_vec{ind_vcv} = 'PEARSON R HAR PER AREA NH H2F'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.P_R2_harvest_2_fish;
 var_string_vec{ind_vcv} = 'PEARSON R2 HAR PER AREA NH H2F';        ind_vcv = ind_vcv + 1;
% var_color_vec{ind_vcv}  = struc.corr.K_harvest_2_fish;
% var_string_vec{ind_vcv} = 'KENDAL TAU HAR PER AREA NH H2F'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_harvest_2_fish;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO HAR PER AREA NH H2F';      ind_vcv = ind_vcv + 1;

 var_color_vec{ind_vcv}  = struc.corr.P_R2_harSA_nh_max5_all;
 var_string_vec{ind_vcv} = 'PEARSON R2 HAR PER AREA SMALL NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
% var_color_vec{ind_vcv}  = struc.corr.K_harSA_nh_max5_all;
% var_string_vec{ind_vcv} = 'KENDAL TAU HAR PER AREA SMALL NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_harSA_nh_max5_all;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO HAR PER AREA SMALL NH MAX5 ALL'; ind_vcv = ind_vcv + 1;

 var_color_vec{ind_vcv}  = struc.corr.P_R2_harMA_nh_max5_all;
 var_string_vec{ind_vcv} = 'PEARSON R2 HAR PER AREA MEDIUM NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
% var_color_vec{ind_vcv}  = struc.corr.K_harMA_nh_max5_all;
% var_string_vec{ind_vcv} = 'KENDAL TAU HAR PER AREA MEDIUM NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_harMA_nh_max5_all;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO HAR PER AREA MEDIUM NH MAX5 ALL'; ind_vcv = ind_vcv + 1;

 var_color_vec{ind_vcv}  = struc.corr.P_R2_harLA_nh_max5_all;
 var_string_vec{ind_vcv} = 'PEARSON R2 HAR PER AREA LARGE NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.corr.K_harLA_nh_max5_all;
%var_string_vec{ind_vcv} = 'KENDAL TAU HAR PER AREA LARGE NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_harLA_nh_max5_all;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO HAR PER AREA LARGE NH MAX5 ALL'; ind_vcv = ind_vcv + 1;

 var_color_vec{ind_vcv}  = struc.corr.P_Sfrac_nh_max5_all;
 var_string_vec{ind_vcv} = 'PEARSON R2 Sfrac NH 1990-2006'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.corr.K_Sfrac_nh_max5_all;
%var_string_vec{ind_vcv} = 'KENDAL TAU Sfrac NH 1990-2006'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_Sfrac_nh_max5_all;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO Sfrac NH 1990-2006'; ind_vcv = ind_vcv + 1;

 var_color_vec{ind_vcv}  = struc.corr.P_Mfrac_nh_max5_all;
 var_string_vec{ind_vcv} = 'PEARSON R2 Mfrac NH 1990-2006'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.corr.K_Mfrac_nh_max5_all;
%var_string_vec{ind_vcv} = 'KENDAL TAU Mfrac NH 1990-2006'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.corr.S_Mfrac_nh_max5_all;
 var_string_vec{ind_vcv} = 'SPEARMAN RHO Mfrac NH 1990-2006'; ind_vcv = ind_vcv + 1;

%-----------------------------------------------------------------------------------------
% bias

%var_color_vec{ind_vcv} = struc.bias.harA_nh_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR PER AREA NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.harA_nh_max_all;
%var_string_vec{ind_vcv} = 'BIAS HAR PER AREA NH MAX ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.harA_top20_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR PER AREA NH TOP20 MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.harSA_nh_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR PER AREA SMALL NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.harMA_nh_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR PER AREA MEDIUM NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.harLA_nh_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR PER AREA LARGE NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.Sfrac_nh_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR Sfrac NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.bias.Mfrac_nh_max5_all;
%var_string_vec{ind_vcv} = 'BIAS HAR Mfrac NH MAX5 ALL'; ind_vcv = ind_vcv + 1;

%-----------------------------------------------------------------------------------------
% RMSE

 var_color_vec{ind_vcv}  = struc.rmse.harA_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR PER AREA NH MAX5 ALL';        ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.rmse.harA_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR PER AREA NH MAX ALL';         ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.rmse.harSA_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR PER AREA SMALL NH MAX5 ALL';  ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.rmse.harMA_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR PER AREA MEDIUM NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.rmse.harLA_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR PER AREA LARGE NH MAX5 ALL';  ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.rmse.Sfrac_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR Sfrac NH MAX5 ALL';           ind_vcv = ind_vcv + 1;
 var_color_vec{ind_vcv}  = struc.rmse.Mfrac_nh_max5_all;
 var_string_vec{ind_vcv} = 'RMSE HAR Mfrac NH MAX5 ALL';           ind_vcv = ind_vcv + 1;

%var_color_vec{ind_vcv} = struc.numvalid.har_S_nohigh_max5_all;
%var_string_vec{ind_vcv} = 'NUMVAL HAR SMALL NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.numvalid.har_M_nohigh_max5_all;
%var_string_vec{ind_vcv} = 'NUMVAL HAR MEDIUM NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.numvalid.har_L_nohigh_max5_all;
%var_string_vec{ind_vcv} = 'NUMVAL HAR LARGE NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.numvalid.Sfrac_nh_max5_all;
%var_string_vec{ind_vcv} = 'NUMVAL HAR Sfrac NH MAX5 ALL'; ind_vcv = ind_vcv + 1;
%var_color_vec{ind_vcv} = struc.numvalid.Mfrac_nh_max5_all;
%var_string_vec{ind_vcv} = 'NUMVAL HAR Mfrac NH MAX5 ALL'; ind_vcv = ind_vcv + 1;

var_color_vec{ind_vcv}  = struc.maxh.year_maxh - 100;
var_string_vec{ind_vcv} = 'YEAR OF MAX HARVEST';                   ind_vcv = ind_vcv + 1;

%-----------------------------------------------------------------------------------------
% list of variables for histograms
%-----------------------------------------------------------------------------------------

ind_vcv_hist = 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.fish;
var_string_vec_hist{ind_vcv_hist} = 'MAX H FISH';           ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.fish_S;
var_string_vec_hist{ind_vcv_hist} = 'MAX H SMALL FISH';     ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.fish_M;
var_string_vec_hist{ind_vcv_hist} = 'MAX H MEDIUM FISH';    ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.fish_L;
var_string_vec_hist{ind_vcv_hist} = 'MAX H LARGE FISH';     ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.harvest;
var_string_vec_hist{ind_vcv_hist} = 'MAX H HARVEST';        ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.harvest_S;
var_string_vec_hist{ind_vcv_hist} = 'MAX H SMALL HARVEST';  ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.harvest_M;
var_string_vec_hist{ind_vcv_hist} = 'MAX H MEDIUM HARVEST'; ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.harvest_L;
var_string_vec_hist{ind_vcv_hist} = 'MAX H LARGE HARVEST';  ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.effort;
var_string_vec_hist{ind_vcv_hist} = 'MAX H EFFORT';         ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.effort_S;
var_string_vec_hist{ind_vcv_hist} = 'MAX H SMALL EFFORT';   ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.effort_M;
var_string_vec_hist{ind_vcv_hist} = 'MAX H MEDIUM EFFORT';  ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist}  = struc.maxh.effort_L;
var_string_vec_hist{ind_vcv_hist} = 'MAX H LARGE EFFORT';   ind_vcv_hist = ind_vcv_hist + 1;

var_color_vec_hist{ind_vcv_hist} = struc.y100.fish;
var_string_vec_hist{ind_vcv_hist} = 'Y100 FISH';            ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.fish_S;
var_string_vec_hist{ind_vcv_hist} = 'Y100 SMALL FISH';      ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.fish_M;
var_string_vec_hist{ind_vcv_hist} = 'Y100 MEDIUM FISH';     ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.fish_L;
var_string_vec_hist{ind_vcv_hist} = 'Y100 LARGE FISH';      ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.harvest;
var_string_vec_hist{ind_vcv_hist} = 'Y100 HARVEST';         ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.harvest_S;
var_string_vec_hist{ind_vcv_hist} = 'Y100 SMALL HARVEST';   ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.harvest_M;
var_string_vec_hist{ind_vcv_hist} = 'Y100 MEDIUM HARVEST';  ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.harvest_L;
var_string_vec_hist{ind_vcv_hist} = 'Y100 LARGE HARVEST';   ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.effort;
var_string_vec_hist{ind_vcv_hist} = 'Y100 EFFORT';          ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.effort_S;
var_string_vec_hist{ind_vcv_hist} = 'Y100 SMALL EFFORT';    ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.effort_M;
var_string_vec_hist{ind_vcv_hist} = 'Y100 MEDIUM EFFORT';   ind_vcv_hist = ind_vcv_hist + 1;
var_color_vec_hist{ind_vcv_hist} = struc.y100.effort_L;
var_string_vec_hist{ind_vcv_hist} = 'Y100 LARGE EFFORT';    ind_vcv_hist = ind_vcv_hist + 1;

nvcv = length(var_color_vec);
nvcvh = length(var_color_vec_hist);

%-----------------------------------------------------------------------------------------
% plot harvest vs. fish in terms of parameters
%-----------------------------------------------------------------------------------------

if (idoprint_scatter_FH)

for indp = 1:nvcv
  
  figs(fig_count) = figure;
  var_color = var_color_vec{indp};
  var_string = var_string_vec(indp);

  scatter(log10(fish),log10(harvest),[5],var_color,'fill');
  cbar_handle = colorbar;
%  set(get(cbar_handle,'ylabel'),'string',var_string,'fontsize',14,'interpreter','none')
  hold on;
  plot(log10(fish),log10(target_harvest_low*ones(size(fish))),'-r','LineWidth',2);
  plot(log10(fish),log10(target_harvest_high*ones(size(fish))),'-r','LineWidth',2);
  hold off;
  ylim([-2 3.5]); xlim([(log10(min(fish))+8) (log10(max(fish))+0.1)]);
  xlabel('log10 FISH AT MAX HARVEST (log10 10^{15} g)','FontSize',14)
  ylabel('log10 MAX HARVEST (log10 10^{12} gwB y^{-1})','FontSize',14)
  title(var_string,'fontsize',14,'interpreter','none')
    
    if findstr(var_string_vec{indp},'PEARSON')
      caxis([0 0.75]);
    end
    if findstr(var_string_vec{indp},'SPEARMAN')
      caxis([0 0.75]);
    end

    if findstr(var_string_vec{indp},'BIAS HAR PER AREA')
      caxis([-5 5]); 
    end
    if findstr(var_string_vec{indp},'BIAS HAR NH')
      caxis([-2e12 2e12]);
    end
    if findstr(var_string_vec{indp},'BIAS HAR SMALL')
      caxis([-2e12 2e12]);
    end
    if findstr(var_string_vec{indp},'BIAS HAR MEDIUM')
      caxis([-2e12 2e12]);
    end
    if findstr(var_string_vec{indp},'BIAS HAR LARGE')
      caxis([-2e12 2e12]);
    end
    if findstr(var_string_vec{indp},'BIAS HAR Sfrac')
      caxis([-0.25 0.25]);
    end
    if findstr(var_string_vec{indp},'BIAS HAR Mfrac')
      caxis([-0.25 0.25]);
    end 
   if findstr(var_string_vec{indp},'RMSE HAR PER AREA')
      caxis([0 5]);
    end
    if findstr(var_string_vec{indp},'RMSE HAR NH')
      caxis([0 5e12]);
    end
    if findstr(var_string_vec{indp},'RMSE HAR SMALL')
      caxis([0 2e12]);
    end
    if findstr(var_string_vec{indp},'RMSE HAR MEDIUM')
      caxis([0 2e12]);
    end
    if findstr(var_string_vec{indp},'RMSE HAR LARGE')
      caxis([0 2e12]);
    end
    if findstr(var_string_vec{indp},'RMSE HAR Sfrac')
      caxis([0 0.5]);
    end
    if findstr(var_string_vec{indp},'RMSE HAR Mfrac')
      caxis([0 0.5]);
    end

  fig_count = fig_count + 1;
end

end 

%-----------------------------------------------------------------------------------------
% scatter plot
% compare satisfactory harvest runs by parameter or variable
% harvest (total or group) vs. variable of interest
%-----------------------------------------------------------------------------------------

if (idoprint_scatter_param_corr)

for indp = 1:1:nvcv
  figs(fig_count) = figure;
  var_color = var_color_vec{indp};
  var_string = var_string_vec(indp);
  scatter(var_color,harvest,[5],'fill')
  ylim([0 500])
  xlabel(var_string,'FontSize',14,'interpreter','none')
  ylabel('MAX HARVEST (10^{12} gwB y^{-1})','FontSize',14)
  % use switch to modify a few of the figures
  hold on
  
  if findstr(var_string_vec{indp},'PEARSON')
    xlim([0 0.75]);  
  end
  if findstr(var_string_vec{indp},'SPEARMAN')
    xlim([0 0.75]);
  end
  
  scatter(var_color(harvest_sat),harvest(harvest_sat),[5],'fill','r')
  hold off

% to make a plot showing satisfactory harvest for all groups and total
% make a Taylor diagram instead (see notes)
%      xlim([0 0.65]); ylim([15 150])
%      h1 = scatter(var_color,harvest,[60],'fill','b');
%      h2 = scatter(var_color(harvest_sat),harvest(harvest_sat),[60],'r','fill');
%      h3 = scatter(var_color(harvest_S_sat),harvest(harvest_S_sat),[30],'g','fill');
%      h4 = scatter(var_color(harvest_M_sat),harvest(harvest_M_sat),[20],'k','fill');
%      h5 = scatter(var_color(harvest_L_sat),harvest(harvest_L_sat),[10],'y','fill');
%      legend([h1 h2 h3 h4 h5],{'ALL','HARVEST SAT','HARVEST S SAT','HARVEST M SAT','HARVEST L SAT'},'Location','northeastoutside');
    
  fig_count = fig_count + 1;
    
end

end

%-----------------------------------------------------------------------------------------
% histograms
%-----------------------------------------------------------------------------------------

if (idoprint_histogram)

for indp = 1:1:nvcv
  figs(fig_count) = figure;
  var_color = var_color_vec{indp};
  var_string = var_string_vec(indp);
  hist(var_color,50)
  ylim([0 500])
%  xlabel(var_string,'FontSize',14,'interpreter','none')
  xlabel(var_string,'FontSize',14)
  ylabel('Frequency','FontSize',14)

  if findstr(var_string_vec{indp},'PEARSON')
    xlim([0 0.75]);
  end
  if findstr(var_string_vec{indp},'SPEARMAN')
    xlim([0 0.75]);
  end
  
% to make a plot showing satisfactory harvest for all groups and total
% make a Taylor diagram instead (see notes)
%      xlim([0 0.65]); ylim([15 150])
%      h1 = scatter(var_color,harvest,[60],'fill','b');
%      h2 = scatter(var_color(harvest_sat),harvest(harvest_sat),[60],'r','fill');
%      h3 = scatter(var_color(harvest_S_sat),harvest(harvest_S_sat),[30],'g','fill');
%      h4 = scatter(var_color(harvest_M_sat),harvest(harvest_M_sat),[20],'k','fill');
%      h5 = scatter(var_color(harvest_L_sat),harvest(harvest_L_sat),[10],'y','fill');
%      legend([h1 h2 h3 h4 h5],{'ALL','HARVEST SAT','HARVEST S SAT','HARVEST M SAT','HARVEST L SAT'},'Location','northeastoutside');
    
  fig_count = fig_count + 1;
    
end  % indp

for indp = 1:1:nvcvh
  figs(fig_count) = figure;
  var_color = var_color_vec_hist{indp};
  var_string = var_string_vec_hist(indp);
  hist(var_color,20)
  ylim([0 500])
  xlabel(var_string,'FontSize',14,'interpreter','none')
  ylabel('Frequency','FontSize',14)

  fig_count = fig_count + 1;
    
end % indp

end

%-----------------------------------------------------------------------------------------
% QQ PLOTS of parameter distributions for satisfactory and unsatisfactory simulations
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Set list of parameters to analyze
param_list = {'E_activation_A','kappa_eppley','Prod_star','te','ppmr','b_allo', ... 
              'E_activation_m','zeta1','h_allo','A00','egg_surv','sel_slope','sel_pos_scale'};

%-----------------------------------------------------------------------------------------
% Set subset of simulations to consider
%subset_s_note = 'MOut.sat.S_A_corr_harvest_all';
subset_s_note = 'MOut.sat.P_A_corr_harvest';
eval(['subset_s = ' subset_s_note ';']);
% unsatisfactory subset of simulations
subset_u = ~subset_s;

if (idoprint_qqplot)

  for indp = 1:1:length(param_list)

    figs(indp) = figure;
    qqplot(MOut.(param_list{indp})(subset_s),MOut.(param_list{indp})(subset_u));
    title([param_list{indp} ' X sat Y unsat'] ,'FontSize',14,'Interpreter','None');

  end % indp
  
end % idoplot_qqplot

%-----------------------------------------------------------------------------------------
% plot time series with color representing a value of the simulation
%-----------------------------------------------------------------------------------------

if (idoprint_timeseries)

% first set of indices
% correlation above 
sat_ind = struc.nruns_vec(struc.sat.harvest);
% array of simulations that satisfy struc.sat.harvest_high
sat_ind_harvest_high = struc.nruns_vec(struc.sat.harvest_high);
% array of simulations that satisfy struc.sat.harvest_low
sat_ind_harvest_low = struc.nruns_vec(struc.sat.harvest_low);

% year_plot_start
year_plot_start = 0;
% year_plot_end
year_plot_end = 300;

% year array
year = struc.year;
Colors = jet;
NoColors = length(Colors);

%-----------------------------------------------------------------------------------------
% time series of biomass and harvest for sat_ind
%-----------------------------------------------------------------------------------------

I = struc.corr.P_R2_harA_nh_max5_all(sat_ind);
Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB = interp1(1:NoColors,Colors,Ireduced);

% biomass time series
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.fish_gi_t(sat_ind,:)*1e-15)
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Biomass (10^{15} gwB)','FontSize',14);
xlim([year_plot_start year_plot_end]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','P R2 HARA','fontsize',14,'interpreter','none')
title('time evolution of biomass for good P R2 HARA runs','fontsize',14);

fig_count = fig_count+1;

% harvest time series
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.harvest_gi_t(sat_ind,:)*struc.parameters.spery*1e-12)
hold on
plot(year,target_harvest_low*ones(size(year)),'-k','LineWidth',3);
plot(year,target_harvest_high*ones(size(year)),'-k','LineWidth',3);
hold off
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Harvest (10^{12} gwB y^{-1})','FontSize',14);
xlim([year_plot_start year_plot_end]); ylim([0 165]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','P R2 HARA','fontsize',14,'interpreter','none')
title('time evolution of harvest for good P R2 HARA runs','fontsize',14);

%-----------------------------------------------------------------------------------------
% time series of biomass and harvest for sat_ind using zeta1 in colorbar
%-----------------------------------------------------------------------------------------

fig_count = fig_count+1;

I = struc.zeta1(sat_ind);
Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB = interp1(1:NoColors,Colors,Ireduced);

% biomass time series
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.fish_gi_t(sat_ind,:)*1e-15)
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Biomass (10^{15} gwB)','FontSize',14);
xlim([year_plot_start year_plot_end]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','zeta1','fontsize',14,'interpreter','none')
title('time evolution of biomass for good P R2 HARA runs','fontsize',14);

fig_count = fig_count+1;

% harvest time series
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.harvest_gi_t(sat_ind,:)*struc.parameters.spery*1e-12)
hold on
plot(year,target_harvest_low*ones(size(year)),'-k','LineWidth',3);
plot(year,target_harvest_high*ones(size(year)),'-k','LineWidth',3);
hold off
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Harvest (10^{12} gwB y^{-1})','FontSize',14);
xlim([year_plot_start year_plot_end]); ylim([0 165]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','zeta1','fontsize',14,'interpreter','none')
title('time evolution of harvest for good P R2 HARA runs','fontsize',14);

%-----------------------------------------------------------------------------------------
% time series of biomass and harvest for sat_ind_harvest_high
%-----------------------------------------------------------------------------------------

fig_count = fig_count+1;

I = struc.corr.P_R2_harA_nh_max5_all(sat_ind_harvest_high);
Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB = interp1(1:NoColors,Colors,Ireduced);

% biomass time series
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.fish_gi_t(sat_ind_harvest_high,:)*1e-15)
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Biomass (10^{15} gwB)','FontSize',14);
xlim([year_plot_start year_plot_end]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','P R2 HARA','fontsize',14,'interpreter','none')
title('time evolution of biomass for high harvest runs','fontsize',14);

fig_count = fig_count+1;

% harvest time series
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.harvest_gi_t(sat_ind_harvest_high,:)*struc.parameters.spery*1e-12)
hold on
plot(year,target_harvest_low*ones(size(year)),'-k','LineWidth',3);
plot(year,target_harvest_high*ones(size(year)),'-k','LineWidth',3);
hold off
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Harvest (10^{12} gwB y^{-1})','FontSize',14);
xlim([year_plot_start year_plot_end]); ylim([0 1050]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','P R2 HARA','fontsize',14,'interpreter','none')
title('time evolution of harvest for high harvest runs','fontsize',14);

%-----------------------------------------------------------------------------------------
% time series of biomass and harvest for sat_ind_harvest_low
%-----------------------------------------------------------------------------------------

fig_count = fig_count+1;

I = struc.corr.P_R2_harA_nh_max5_all(sat_ind_harvest_low);
Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB = interp1(1:NoColors,Colors,Ireduced);

% plot fish biomass time series for sat.harvest_low runs
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.fish_gi_t(sat_ind_harvest_low,:)*1e-15)
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Biomass (10^{15} gwB)','FontSize',14);
xlim([year_plot_start year_plot_end]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','P R2 HARA','fontsize',14,'interpreter','none')
title('time evolution of biomass for low harvest runs','fontsize',14);

fig_count = fig_count+1;

% plot harvest time series for sat.harvest_low runs
figs(fig_count) = figure;
set(gcf,'DefaultAxesColorOrder',RGB)
plot(year,struc.harvest_gi_t(sat_ind_harvest_low,:)*struc.parameters.spery*1e-12)
hold on
plot(year,target_harvest_low*ones(size(year)),'-k','LineWidth',3);
plot(year,target_harvest_high*ones(size(year)),'-k','LineWidth',3);
hold off
colorbar
caxis([min(I),max(I)])
xlabel('time','FontSize',14); ylabel('Harvest (10^{12} gwB y^{-1})','FontSize',14);
xlim([year_plot_start year_plot_end]); ylim([0 100]);
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'string','P R2 HARA','fontsize',14,'interpreter','none')
title('time evolution of harvest for low harvest runs','fontsize',14);

end

%-----------------------------------------------------------------------------------------
% print figures to file
%-----------------------------------------------------------------------------------------

if (idoprint)

  if (idoprint_scatter_FH)
    filename = 'fig_MCV3_4_FH_scatter.ps';
  elseif (idoprint_scatter_param_corr)
    filename = 'fig_MCV3_4_param_corr_scatter.ps';
  elseif (idoprint_histogram)
    filename = 'fig_MCV3_4_histogram.ps';
  elseif (idoprint_timeseries)
    filename = 'fig_MCV3_4_timeseries.ps';
  elseif (idoprint_qqplot)
    filename = 'fig_MCV3_4_analyze_qqplot.ps';
  end

  if (idoprint_scatter_FH | idoprint_scatter_param_corr | idoprint_histogram | idoprint_timeseries | idoprint_qqplot)

    for idx = 1:length(figs)
      print(figs(idx), '-append', '-dpsc2',filename);
    end

  end % OR statements

end % idoprint

%-----------------------------------------------------------------------------------------
% END OF SCRIPT
