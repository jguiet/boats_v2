function [boats] = run_boats(Ecological_frc, Economical_frc, Userdef_Params)
%**************************************************************************
% BOATS MAIN
% The BiOeconomic mArine Trophic Size-spectrum model computes the
% harvestable biomass of fish. Coupled with a bioeconomical module it
% simulates the dynamic of fishing. The governing equations of the model
% are described in:
% - Carozza, D. A.; Bianchi, D. & Galbraith, E. D. 
% The ecological module of BOATS-1.0: a bioenergetically-constrained model 
% of marine upper trophic levels suitable for studies of fisheries and 
% ocean biogeochemistry. 
% Geoscientific Model Development, 2015, 8, 10145-10197 
% - Carozza, D. A.; Bianchi, D. & Galbraith, E. D. 
% Formulation, General Features and Global Calibration of a 
% Bioenergetically-Constrained Fishery Model.
% PLOS ONE, Public Library of Science, 2017, 12, 1-28
%**************************************************************************
% Inputs:
%   - parameters.txt
%   - Ecological.mat
%   - Economical.mat (for harvest simulations)
% Outputs:
%   - Boats.mat
%**************************************************************************
 addpath('main')
 addpath('general')
 
 
%**************************************************************************
% LOAD INPUT
%**************************************************************************
 parameters;
 
 
%**************************************************************************
% PREPARE SIMULATION
%**************************************************************************
% Forcings datasets *******************************
 forcing_ecological=Ecological_frc;
 forcing_economical=Economical_frc;
% Make output/restar dirs *************************
 if ~exist(boats.param.path.outdir)
    mkdir(boats.param.path.outdir)
 end

 
%**************************************************************************
% IF ENSEMBLE PARAMETERS
%**************************************************************************
 if boats.param.main.param_ens 
    load(Userdef_Params);
    param_loops=length(ens_index);
 else
    param_loops=1;
 end

 
%**************************************************************************
% MAIN LOOP on parameter sets
%**************************************************************************
 for indr = 1:param_loops     
   %***********************************************************************
   % UPDATE PARAMETERS FOR ENSEMBLES
   %***********************************************************************
   if boats.param.main.param_ens
     boats.param.main.sname_rest = ['_ind_' num2str(ens_index(indr))];
     boats.param = modify_parameters(boats,...
                    'E_activation_A',ens_param.E_activation_A(indr,:), ...
                    'kappa_eppley',ens_param.kappa_eppley(indr,:), ...
                    'Prod_star',ens_param.Prod_star(indr,:), ...
                    'te',ens_param.te(indr,:),...
                    'ppmr',ens_param.ppmr(indr,:), ...
                    'tro_sca',ens_param.tro_sca(indr,:), ...
                    'b_allo', ens_param.b_allo(indr,:), ...
                    'E_activation_m',ens_param.E_activation_m(indr,:), ...
                    'zeta1',ens_param.zeta1(indr,:),...
                    'h_allo',ens_param.h_allo(indr,:), ...
                    'A00',ens_param.A00(indr,:),...
                    'egg_surv',ens_param.egg_surv(indr,:), ...
                    'sel_slope',ens_param.sel_slope(indr,:), ...
                    'sel_pos_scale',ens_param.sel_pos_scale(indr,:));
   else
     boats.param.main.sname_rest = [''];
   end

   
   %***********************************************************************
   % LOAD_FORCINGS
   %***********************************************************************   
   boats.forcing = load_forcing(boats,...
                     forcing_ecological,forcing_economical);
   
                 
   %***********************************************************************
   % PREPARE OUTPUT
   %***********************************************************************  
   boats.output = initialize_output({'annual'});
   
   
   %***********************************************************************
   % PREPARE STRUCTURE
   %***********************************************************************  
   boats.structure = set_structure(boats);           
   
   
   %***********************************************************************
   % INITIALIZE SPECTRUM
   %***********************************************************************  
   boats.initial = initialize_domains(boats);

   
   %***********************************************************************
   % INTEGRATE BOATS
   %***********************************************************************
   disp(['Running BOATS # ' boats.param.main.sim_name '_' boats.param.main.sim_type boats.param.main.sname_rest]); 
   boats = integrate(boats);                                  

   %-----------------------------------------------------------------------------------------------------------
   % Remove unnecessary files and save restart
   %-----------------------------------------------------------------------------------------------------------
   if boats.param.main.save_restart==1
     cd(boats.param.path.outdir)
     boats = save_restart(boats);
     cd(boats.param.path.wrkdir)
   end

   %------------ -----------------------------------------------------------------------------------------------
   % Save output
   %-----------------------------------------------------------------------------------------------------------
   if boats.param.main.save_output==1
     savenameall = [boats.param.main.sim_name '_' boats.param.main.sim_type boats.param.main.sname_rest '.mat'];
     disp(['saving boats2d structure:' savenameall]);
     parsave_boats([boats.param.path.outdir savenameall],boats);                            % Note: parsave is for parallel runs   
   end
 end
%--------------------------------------------------------------------------------------------------------------
% END OF MAIN LOOP
%--------------------------------------------------------------------------------------------------------------


%**************************************************************************************************************
% END OF MAIN SCRIPT
