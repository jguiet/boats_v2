% Plot the timeseries for the best ensembles
clear all

% Load data
load('listbest.mat')

% Plot
for k = 1:50
	name = strcat('/data/results/jguiet/BOATS/MCRun/data/',list(k))
	load(name{1})
	B(k,:) = boats.output.annual.fish_gi_t';
	H(k,:) = boats.output.annual.harvest_gi_t';
	E(k,:) = boats.output.annual.effort_gi_t';
end

% Save
save('BHE.mat','B','H','E')


