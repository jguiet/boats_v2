% analyze_boats2d_MC_sat_unsat_stats.m
%-----------------------------------------------------------------------------------------
% compares basic statistical properties of the random variates (parameters)
% of the best runs of a Monte-Carlo simulation to the complete set
% in order to determine whether the distribution is different
%-----------------------------------------------------------------------------------------
% kstest2
% Two-sample Kolmogorov-Smirnov test with null hypothesis that both samples are from the
% same distribution. Result is 1 if null hypothesis is rejected.

%-----------------------------------------------------------------------------------------
% load processed structure of Monte Carlo simulations
 load /archive/dcarozza/boats2d_MC/processed/MOut_MCV3_5.mat

%-----------------------------------------------------------------------------------------
% load parameter values
 load /archive/dcarozza/boats2d_MC/processed/parameters_mc.mat
 [nparams nruns] = size(parameters_mc);

%-----------------------------------------------------------------------------------------
% add column for tau (log10(te) / log10(ppmr)
parameters_mc(14,:) = log10(parameters_mc(4,:)) ./ log10(parameters_mc(5,:));

%-----------------------------------------------------------------------------------------
% calculate satisfactory and unsatisfactory subsets of simulations in order to calculate
% test and goodness of fit statistics

%-----------------------------------------------------------------------------------------
% satisfactory subset of simulations
subset_s_note = 'total harvest between 70 and 150, no h2f ratio above 0.4, medium H between 0.3 and 1.2 small, large H between 0.1 and 0.8 S, Pearson R2 above 0.45';

subset_s = (MOut.maxh.harvest < 150) & (MOut.maxh.harvest > 70) & ...
            (MOut.h2f.max < 0.4) & ...
            (MOut.maxh.harvest_M > 0.3*MOut.maxh.harvest_S) & (MOut.maxh.harvest_M < 1.2*MOut.maxh.harvest_S)  & ...
            (MOut.maxh.harvest_L > 0.1*MOut.maxh.harvest_S) & (MOut.maxh.harvest_L < 0.8*MOut.maxh.harvest_S) & ...
            (MOut.corr2006.P_R2_harA_nh_max_all > 0.45);

% unsatisfactory subset of simulations
subset_u = ~subset_s;

%-----------------------------------------------------------------------------------------
% Set list of parameters to analyze
param_list = {'E_activation_A','kappa_eppley','Prod_star','te','ppmr','b_allo', ... 
              'E_activation_m','zeta1','h_allo','A00','egg_surv','sel_slope','sel_pos_scale','tau'};

%-----------------------------------------------------------------------------------------
% Calculate mean, median, std
% Calculate test and goodness of fit statistics

for indp = 1:1:(nparams+1)
  
  % mean of satisfactory subset
  struc.(param_list{indp}).mean_s       = mean(parameters_mc(indp,subset_s));
  % mean of unsatisfactory subset
  struc.(param_list{indp}).mean_u       = mean(parameters_mc(indp,subset_u));

  % median of satisfactory subset
  struc.(param_list{indp}).median_s     = median(parameters_mc(indp,subset_s));
  % median of unsatisfactory subset
  struc.(param_list{indp}).median_u     = median(parameters_mc(indp,subset_u));

  % standard deviation of satisfactory subset
  struc.(param_list{indp}).std_s        = std(parameters_mc(indp,subset_s));
  % standard deviation of unsatisfactory subset
  struc.(param_list{indp}).std_u        = std(parameters_mc(indp,subset_u));

  % rejection (1) or acceptance (0) of null hypothesis of 2-sample t-test (comparing means)
  struc.(param_list{indp}).ttest2_s_u   = ttest2(parameters_mc(indp,subset_s),parameters_mc(indp,subset_u));
  % rejection (1) or acceptance (0) of null hypothesis of 2-sample F-test (comparing variances)
  struc.(param_list{indp}).vartest2_s_u = vartest2(parameters_mc(indp,subset_s),parameters_mc(indp,subset_u));
  % rejection (1) or acceptance (0) of null hypothesis of 2-sample Kolmogorov-Smirnov test (comparing distributions)
  % p-value (2-sided) of 2-sample Kolmogorov-Smirnov test (comparing distributions)
  [struc.(param_list{indp}).kstest2_s_u_p05 struc.(param_list{indp}).kstest2_s_u_pval] = kstest2(parameters_mc(indp,subset_s),parameters_mc(indp,subset_u));
  
end

% include additional variables to results structure
%struc.subset_s_note = subset_s_note;
struc.subset_s = subset_s;
struc.subset_u = subset_u;

% rename results structure
analysis_stats_boats2d_MC = struc;

% save result structure
save /archive/dcarozza/writing/01_in_preparation/cbg_boats_econ_plosone/output/analysis_stats_boats2d_MC.mat analysis_stats_boats2d_MC -v7.3

%-----------------------------------------------------------------------------------------
% end of script