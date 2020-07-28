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

% Loop on parameter sets
ens_param = [];
for indr = 1:param_loops
    if (boats.param.ecology.demersal)&&(boats.param.ecology.pelagic)
        ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic' , indr,ens_param);
        ens_param=load_MC(boats.param.main.MC_param_Dem,'Demersal', indr,ens_param);
    else
        ens_param=load_MC(boats.param.main.MC_param_Pel,'Pelagic',[]);
    end
end 
save('ens_param_A.mat','ens_param')

% Launch BOATS in serial on nrun replicates
run_boats('frc/Ecological.mat','frc/Economical.mat',ens_param)
