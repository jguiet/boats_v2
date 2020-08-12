%**************************************************************************************************************
% FUNCTION set_MC.m
% Determine Monte Carlo random parameters ranges for the simulation
%**************************************************************************************************************
function MParams = set_MC(Community)

%--------------------------------------------------------
% Note: need to use rand to generate random numbers
% and be able to re-set a different seed each run
% Re-seeds the random number generator stresm
 rstate = round(sum(clock));
 disp(['Random Seed: ' num2str(rstate)]);
 rand('state',rstate);


%--------------------------------------------------------
 % Montecarlo simulation parameters
 % structure is: 
 % {'name', 'mean', 'std', 'distribution'}

 % zeta1
 % Gislason et al. (2010) range for parameter 'a', need to transform this to exp(a)
 % log normally distributed
 zeta1_mean = 0.55;
 zeta1_std  = 1.12 - zeta1_mean;
 
 % h_allo
 % Gislason et al. (2010) range
 % Quantities in paper are for length (mass is to the 1/3 and so need to divide by 3
 h_allo_mean = 1.61/3;
 h_allo_std  = (1.61-1.34)/3;
  
 % egg_surv
 % from 10^(-3.5) to 0.05
 egg_surv_low  = 10^(-3.5);
 egg_surv_high = 0.10;
 egg_surv_mean = (egg_surv_low + egg_surv_high)/2;
 egg_surv_std  = (egg_surv_high - egg_surv_mean)/sqrt(3);

 % sel_slope
 % uniform from 12 to 24
 sel_slope_low  = 12;
 sel_slope_high = 24;
 sel_slope_mean = (sel_slope_low + sel_slope_high)/2;
 sel_slope_std  = (sel_slope_high - sel_slope_mean)/sqrt(3);

 % sel_pos_scale
 % uniform from 0.5 to 1.5
 sel_pos_scale_low  = 0.5;
 sel_pos_scale_high = 1.5;
 sel_pos_scale_mean = (sel_pos_scale_low + sel_pos_scale_high)/2;
 sel_pos_scale_std  = (sel_pos_scale_high - sel_pos_scale_mean)/sqrt(3); 

 % mc_benthic
 mc_benthic_low = 1/5000*0.001;
 mc_benthic_high = 1/5000*0.5;
 mc_benthic_mean = (mc_benthic_low + mc_benthic_high)/2;
 mc_benthic_std  = (mc_benthic_high - mc_benthic_mean)/sqrt(3);

 if strcmp(Community,'Pelagic')
 MParams = {
		'E_activation_A',	 0.45,	             0.09,          NaN,   NaN,     NaN,     'normal';
        'E_activation_m',    0.45,	             0.09,          NaN,   NaN,     NaN,     'normal';
        'b_allo',		     0.7,	             0.2,           NaN,   NaN,     NaN,     'normal';
        'A00',               4.46,	             0.5,           NaN,   NaN,     NaN,     'normal';
        'te',                0.13,	             0.04,          NaN,   NaN,     NaN,     'uniform';
	    'ppmr',              5000,	             2500,          NaN,   NaN,     NaN,     'uniform';
%	    'kappa_eppley',	     0.0631,             0.009,	        NaN,   NaN,     NaN,     'normal';
%		'Prod_star',	     0.37,	             0.1,	        NaN,   NaN,     NaN,     'normal';        
        'zeta1',   	         zeta1_mean,	     zeta1_std,     NaN,   NaN,     NaN,     'normal';
%        'zeta1',	     0.45,	 0.24,	 'lognormal';
        'h_allo',	         h_allo_mean,        h_allo_std,    NaN,   NaN,     NaN,     'normal';
%        'h_allo',	  0.44,		0.18,		'normal';
        'egg_surv',          0.05,               0.028,         NaN,   NaN,     NaN,     'uniform';
%        'sel_slope',         sel_slope_mean,     sel_slope_std, NaN,   NaN,     NaN,     'uniform';
        'sel_pos_scale',     0.8,                0.288,         NaN,   NaN,     NaN,     'uniform';
        'mc_benthic',        mc_benthic_mean,    mc_benthic_std,NaN,   NaN,     NaN,     'uniform'};  
 elseif strcmp(Community,'Demersal')
 MParams = {
		'E_activation_A',	 0.45,	             0.09,          NaN,   NaN,     NaN,     'normal';
        'E_activation_m',    0.45,	             0.09,          NaN,   NaN,     NaN,     'normal';
        'b_allo',		     0.7,	             0.2,           NaN,   NaN,     NaN,     'normal';
        'A00',               4.46,	             0.5,           NaN,   NaN,     NaN,     'normal';
        'te',                0.13,	             0.04,          NaN,   NaN,     NaN,     'uniform';
	    'ppmr',              5000,	             2500,          NaN,   NaN,     NaN,     'uniform';
%	    'kappa_eppley',	     0.0631,             0.009,	        NaN,   NaN,     NaN,     'normal';
%		'Prod_star',	     0.37,	             0.1,	        NaN,   NaN,     NaN,     'normal';        
        'zeta1',   	         zeta1_mean,	     zeta1_std,     NaN,   NaN,     NaN,     'normal';
%        'zeta1',	     0.45,	 0.24,	 'lognormal';
        'h_allo',	         h_allo_mean,        h_allo_std,    NaN,   NaN,     NaN,     'normal';
%        'h_allo',	  0.44,		0.18,		'normal';
        'egg_surv',          0.05,               0.028,         NaN,   NaN,     NaN,     'uniform';
%        'sel_slope',         sel_slope_mean,     sel_slope_std, NaN,   NaN,     NaN,     'uniform';
        'sel_pos_scale',     0.8,                0.288,         NaN,   NaN,     NaN,     'uniform';  
        'mc_benthic',        mc_benthic_mean,    mc_benthic_std,NaN,   NaN,     NaN,     'uniform'};
 end

%**************************************************************************************************************
% END FUNCTION

