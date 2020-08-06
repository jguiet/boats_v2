% **********************************************
% This script launches BOATS in Monte Carlo mode
% **********************************************
clear all

addpath('mc/');
parameters

% Set parameter distributions for undetermined distributions
boats.param.main.MC_param_Pel=set_MC('Pelagic');
boats.param.main.MC_param_Dem=set_MC('Demersal');
param_loops=boats.param.main.nrun;

% Loop ensemble of simulations -> spin-up + transient
for indr = 1:param_loops
    ens_param = [];
    % Parameter selection
    if (boats.param.ecology.demersal)&&(boats.param.ecology.pelagic)
        ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic' , 1,ens_param);
        ens_param=load_MC(boats.param.main.MC_param_Dem,'Demersal', 1,ens_param);
    else
        ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic',[]);
    end
    % run
    run_boats('Ecological.mat','Economical_q7.mat',ens_param,'PP','nh')
    run_boats('Ecological.mat','Economical_q7.mat',ens_param,'restart','h')
end 

