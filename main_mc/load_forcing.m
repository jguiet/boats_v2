%**************************************************************************************************************
% FUNCTION load_forcing.m
% Load forcings preprocessed by "preprocess.m" for the simulation :
% Ecological.mat
% Economical.mat
%**************************************************************************************************************
function forcing = load_forcing(boats,forcing_ecological,forcing_economical)

%---------------------------------
% Forcing ecological:
if exist(forcing_ecological,'file')
    load(forcing_ecological);
else
    disp('Hum, double-check the path for ecological forcing:')
    disp(forcing_ecological)
end
forcing.mask=repmat(Ecological.mask,[1 1 size(Ecological.npp,3)]);
forcing.nlat=size(forcing.mask,1);
forcing.nlon=size(forcing.mask,2);
forcing.npp=Ecological.npp;
forcing.npp(find(forcing.mask==1))=NaN;
forcing.npp_ed=Ecological.npp_ed;
forcing.npp_ed(find(forcing.mask==1))=NaN;
forcing.pfb=Ecological.pfb;
forcing.pfb(find(forcing.mask==1))=NaN;
forcing.temperature_pel=Ecological.temperature_pel;
forcing.temperature_pel_K=Ecological.temperature_pel+boats.param.conversion.C_2_K;
forcing.temperature_dem=Ecological.temperature_dem;
forcing.temperature_dem_K=Ecological.temperature_dem+boats.param.conversion.C_2_K;
forcing.surf=Ecological.surface;

%--------------------------------- 
% Forcing economical
if strcmp(boats.param.main.sim_type,'h')
    if exist(forcing_economical,'file')
        load(forcing_economical);
    else
        disp('Hum, double-check the path for economical forcing:')
        disp(forcing_economical)
    end
    load(forcing_economical)
    forcing.price=Economical.price;
    forcing.cost=Economical.cost;
    forcing.catchability=Economical.catchability;
end
 
%---------------------------------
  % Convert maps to vectors
  for itime = 1:size(forcing.npp,3)
      [forcing.npp_vec(:,itime) forcing.indlat forcing.indlon]           = function_map_2_vec(squeeze(forcing.npp(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.npp_ed_vec(:,itime) forcing.indlat forcing.indlon]        = function_map_2_vec(squeeze(forcing.npp_ed(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.pfb_vec(:,itime) forcing.indlat forcing.indlon]           = function_map_2_vec(squeeze(forcing.pfb(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_pel_vec(:,itime) forcing.indlat forcing.indlon]   = function_map_2_vec(squeeze(forcing.temperature_pel(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_pel_K_vec(:,itime) forcing.indlat forcing.indlon] = function_map_2_vec(squeeze(forcing.temperature_pel_K(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_dem_vec(:,itime) forcing.indlat forcing.indlon]   = function_map_2_vec(squeeze(forcing.temperature_dem(:,:,itime)),squeeze(forcing.mask(:,:,1)));
      [forcing.temperature_dem_K_vec(:,itime) forcing.indlat forcing.indlon] = function_map_2_vec(squeeze(forcing.temperature_dem_K(:,:,itime)),squeeze(forcing.mask(:,:,1)));
  end % itime
  [forcing.surf_vec forcing.indlat forcing.indlon]                = function_map_2_vec(forcing.surf,squeeze(forcing.mask(:,:,1)));
  forcing.nvec=size(forcing.surf_vec,1);
%**************************************************************************************************************
% END FUNCTION

