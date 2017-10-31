 function boats = save_restart(boats)

 % Saves the required restart
 filename = boats.param.main.sname_rest;
 restart.dfish 	= squeeze(boats.dfish);
 dir_restart = boats.param.path.outdir;

 if strcmp(boats.param.main.sim_type,'h')
    restart.effort = squeeze(boats.effort);
 end
 
 % Here puts other things that can be useful for restart - e.g. run parameters
 if (0)
    restart.parameters = boats.parameters;
 end

 %-------------------------------------------------------------------------------
 % remove unneeded arrays
 %-------------------------------------------------------------------------------
 boats.parameters.data_monthly = [];
 
 if isfield(boats,'dfish')
   disp(['remove dfish from BOATS']);
   boats = rmfield(boats,'dfish');
 end
 if isfield(boats,'fish_gi_g_temp')
   disp(['remove fish_gi_g_temp from BOATS']);
   boats = rmfield(boats,'fish_gi_g_temp');
 end
 if isfield(boats,'fish_sel_gi_temp')
   disp(['remove fish_sel_gi_temp from BOATS']);
   boats = rmfield(boats,'fish_sel_gi_temp'); 
 end
 if isfield(boats,'ESM_forcing')
   disp(['remove ESM_forcing from BOATS']);
   boats = rmfield(boats,'ESM_forcing');
 end
 
 if strcmp(boats.param.main.sim_type,'h')
   
   if isfield(boats,'effort')
     disp(['remove effort from BOATS']);
     boats = rmfield(boats,'effort');
   end
   if isfield(boats,'effort_gi_g_temp')
     disp(['remove effort_gi_g_temp from BOATS']); 
     boats = rmfield(boats,'effort_gi_g_temp'); 
   end
   if isfield(boats,'harvest_gi_g_temp')
     disp(['remove harvest_gi_g_temp from BOATS']);
     boats = rmfield(boats,'harvest_gi_g_temp'); 
   end
   if isfield(boats,'fish_gi_temp')
     disp(['remove fish_gi_temp from BOATS']); 
     boats = rmfield(boats,'fish_gi_temp'); 
   end
   if isfield(boats,'fish_sel_gi_temp')
     disp(['remove fish_sel_gi_temp from BOATS']); 
     boats = rmfield(boats,'fish_sel_gi_temp'); 
   end
   if isfield(boats,'harvest_gi_temp')
     disp(['remove harvest_gi_temp from BOATS']); 
     boats = rmfield(boats,'harvest_gi_temp'); 
   end
   if isfield(boats,'effort_gi_temp')
     disp(['remove effort_gi_temp from BOATS']); 
     boats = rmfield(boats,'effort_gi_temp'); 
   end
   
 end

 % display file to save
 disp(['saving restart file: ' dir_restart '/' 'restart_' boats.param.main.sim_name '_' boats.param.main.sim_type  boats.param.main.sname_rest '.mat']);
 % save boats restart
 save([dir_restart '/' 'restart_' boats.param.main.sim_name '_' boats.param.main.sim_type filename],'restart');

 % add restart to the boats structure
 %boats.restart = restart;