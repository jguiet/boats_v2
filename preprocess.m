%**************************************************************************
% BOATS PRE-PROCESS
% Definition of the forcing for the BOATS model ('run_boats.m'). 
% This preprocessing step allows to load and convert different forcing
% sources (Observations/ROMS/...) in order to convert them into suitable
% inputs for BOATS.
%**************************************************************************
% BOATS inputs generated :
% Ecological input 'Ecological.mat' (for no harvest & harvest simulations):
%   - mask [boolean] (mask of land = 1 and ocean = 0 cells in the domain)
%   - lon/lat [degrees] (coordinates of the computation nodes)
%   - surface [m^2] (surface of the computation cells for integration)
%   - npp [mmolC m^-2 s^-1] (net primary production)
%   - npp_ed [mmolC m^-3 d^-1] (net primary production averaged on euphotic zone depth)
%   - temperature [deg C] (water temperature)
% Economical input 'Economical.mat' (for harvest simulations):
%   - price [$ g^-1] (fish price history)
%   - cost [$ W^-1] (exploitation cost history)
%   - catchability [m^2 W^-1 s^-1] (catchability history)
%   - efftar [s-1] (effective effort target)
%   - societenf [??] (societal enforcement)
%**************************************************************************
addpath('preprocess')
clear all


%**************************************************************************
% DEFINE FORCING CHARACTERISTICS
%**************************************************************************
% Preprocess options *******************************
plot_input = 1;                                     % (yes 1 or no 0)
create_ecology = 1;                                 % (yes 1 or no 0)
create_economy = 1;                                 % (yes 1 or no 0)
create_regulation = 1;                              % (yes 1 or no 0)

% General forcing paths and characteristics ********
nlat = 240;                                         % nlat = m nodes along latitude
nlon = 120;                                         % nlon = n nodes along longitude
ntime = 12;                                         % ntime = t distinct time steps per forcing
ngroup = 3;                                         % ngroup = g fish groups per forcing (facultativ)
nensemble = 5;                                      % nensemble = e ensembles per forcing (facultativ)
% mask
mask_path = 'input/BioForcing/ROMS181v17MonthlyMeans.mat'; % Path of forcing dataset where is the mask
mask_var = 'ACmodelMonthlyMeans.landmask';          % Name of mask variable in mask_path
wet=1;                                              % Value for ocean cells in mask_var
dry=NaN;                                            % Value for continent cells in mask_var
mask_dim = [nlat nlon 1 1 1];                       % Dimension of the mask forcing generated
% coordinates lon/lat
lon_path = 'input/BioForcing/ROMS181v17MonthlyMeans.mat';  % Path of forcing dataset where is the longitude
lon_var = 'ACmodelMonthlyMeans.lonr';               % Name of longitude variable in lon_var
lat_path = 'input/BioForcing/ROMS181v17MonthlyMeans.mat';  % Path of forcing dataset where is the latitude
lat_var = 'ACmodelMonthlyMeans.latr';               % Name of latitude variable in lat_var
lon_dim = [nlat nlon 1 1 1];                        % Dimension of the longitude array generated
lat_dim = [nlat nlon 1 1 1];                        % Dimension of the latitude array generated
% surface
surf_path = 'input/BioForcing/ACM075_120_240_area.mat';    % Path of forcing dataset where is the surface
surf_var = 'ACM075_120_240_area';                   % Name of surface variable in surf_path
surf_unit = '[km^2]';                               % surf_var unit ([km^2] or [m^2] or ??)
surf_dim = [nlat nlon 1 1 1];                       % Dimension of the surface array generated

% Ecological forcing paths and characteristics *****
% depth
depth_path = 'input/BioForcing/ROMS181v17MonthlyMeans.mat';% Path of forcing dataset where is the depth
depth_var = 'ACmodelMonthlyMeans.bathy';            % Name of depth variable in depth_path
depth_type = 'value';                               % Type of depth definition (map or value)
depth_unit = '[m]';                                 % depth_var unit ([m] or [km] or ??)
ed=75;                                              % User defined depth of euphotic zone
depth_dim = [nlat nlon 1 1 1];                      % Dimension of the depth array generated
% npp
npp_path = 'input/BioForcing/ROMS181v17MonthlyMeans.mat';  % Path of forcing dataset where is the primary production
npp_var = 'ACmodelMonthlyMeans.ppC';                % Name of primary porduction variable in npp_path
npp_unit = '[mgC m^-2 d^-1]';                       % npp_var unit ([mgC m^-2 d^-1] or ??)
npp_dim = [nlat nlon ntime 1 1];                    % Dimension of the npp array generated
% temperature
temp_path = 'input/BioForcing/ROMS181v17MonthlyMeans.mat'; % Path of forcing dataset where is the temperature
temp_var = 'ACmodelMonthlyMeans.bottomT';           % Name of temperature variable in temp_path         
temp_unit = '[degC]';                               % temp_var unit ([degC])
temp_dim = [nlat nlon ntime 1 1];                   % Dimension of the temp array generated

% Economical forcing paths and characteristics *****
% price
price_path = 'input/EcoForcing/price_forcing.mat'; % Path of forcing dataset where is the price
price_var = 'udef';                                 % Name of price variable in price_path or user defined (udef)
price_type = 'cst';                                 % Type of user defined price if price_var = 'udef' (constant cst or ??)
price_ref1 = 1.264079636832035e-04;                 % First parameter for user defined price forcing
price_dim=[1 1 3600];                               % Dimension of user defined price forcing
price_unit = '[$ g^-1]';                            % price_var unit ([$ g^-1] or ??)
% cost
cost_path = 'input/EcoForcing/cost_effort_forcing.mat'; % Path of forcing dataset where is the cost
cost_var = 'udef';                                  % Name of cost variable in cost_path or user defined (udef)
cost_type = 'cst';                                  % Type of user defined cost if cost_var = 'udef' (constant cst or ??)
cost_ref1 = 1.8520e-07;                             % First parameter for user defined cost forcing
cost_dim=[1 1 3600];                                % Dimension of user defined cost forcing
cost_unit = '[$ W^-1]';                             % cost_var unit ([$ W^-1] or ??)
% catchability
catch_path = 'input/EcoForcing/catchability_forcing.mat'; % Path of forcing dataset where is the catchability
catch_var = 'udef';                                 % Name of catchability variable in catch_path or user defined (udef)
catch_type = 'rate';                                % Type of user defined cost if cost_var = 'udef' (constant cst or rate or ??)
catch_ref1 = 7.6045e-08;                            % First parameter for user defined catchability forcing
catch_ref2 = 0.05;                                  % Second parameter for user defined catchability forcing
catch_dim=[1 1 3600];                               % Dimension of user defined catchability forcing
catch_unit = '[m^2 W^-1 s^-1]';                     % catch_var unit ([m^2 W^-1 s^-1] or ??)

% Regulation forcing paths and characteristics *****
% MSY effective effort target
efftarg_path = '/Users/jguiet/MODELS/BOATS/forcing/dev/E_eff_MSY_true_all.mat'; % Path of forcing dataset where is t??
efftarg_var = 'E_eff_MSY_true_all';
efftarg_dim = [nlat nlon 1 ngroup nensemble];       % Dimension of the effort target array generated
efftarg_unit = '[s-1 (?)]';                         % efftarg_var unit ([s^-1])
% Social enforcement target
% GOOD LUCK...(MODIF)
%**************************************************************************
% END DEFINE FORCING CHARACTERISTICS
%**************************************************************************


%**************************************************************************
% LOAD and CONVERT
%**************************************************************************
disp('Lets go !!!')

% Load and convert ecological forcing *************
if create_ecology
    disp('Create ecological forcing')
    % Create mask *********************************
    mask=get_var(mask_path,mask_var,mask_dim);
    mask=arrange_mask(mask,wet,dry);
    % Create coordinates lon/lat ******************
    lon=get_var(lon_path,lon_var,lon_dim);
    lat=get_var(lat_path,lat_var,lat_dim);
    %ATT JG [lon lat]=arrange_coord(lon,lat);
    % Create surface ******************************
    surface=get_var(surf_path,surf_var,surf_dim);
    surface=arrange_surf(surface,surf_unit,lon,lat);
    % Create npp *********************************
    npp=get_var(npp_path,npp_var,npp_dim);
    depth=get_var(depth_path,depth_var,depth_dim);
    [npp npp_ed]=arrange_npp(npp,npp_unit,depth,depth_unit,npp_dim,ed,depth_type);
    % Create temperature *************************
    temperature=get_var(temp_path,temp_var,temp_dim);
    % Plot forcings ******************************
    if plot_input
        % surface
        plot_domain2D(surface,lon,lat,mask,'surface [m^2]',1)
        % npp
        plot_domain2D(npp,lon,lat,mask,'npp [mmolC m^{-2} s^{-1}]',2)
        plot_domain2D(npp_ed,lon,lat,mask,'npp_{ed} [mmolC m^{-3} d^{-1}]',4)
        % temperature
        plot_domain2D(temperature,lon,lat,mask,'temperature [deg C]',3)
    end
    % Save forcing
    Ecological.mask=mask;
    Ecological.lon=lon;
    Ecological.lat=lat;
    Ecological.surface=surface;
    Ecological.npp=npp;
    Ecological.npp_ed=npp_ed;
    Ecological.temperature=temperature;
    save('Ecological.mat','Ecological','-v7.3')
end

% Load and convert economical forcing *************  
if create_economy
    disp('Create economical forcing')
    % Create price ********************************
    price=get_var(price_path,price_var,[]);
    price=udef_var(price_type, price_dim, price_ref1);
    % ATT JG arrange
    % Create price ********************************
    cost=get_var(cost_path,cost_var,[]);
    cost=udef_var(cost_type, cost_dim, cost_ref1);
    % ATT JG arrange
    % Create price ********************************
    catchability=get_var(catch_path,catch_var,[]);
    catchability=udef_var(catch_type, catch_dim, catch_ref1, catch_ref2);
    % ATT JG arrange 
    % Plot forcings ******************************
    if plot_input
        % price
        plot_domain1D(1:length(price),price,'time','price [$ g^-1]',5)
        % cost
        plot_domain1D(1:length(cost),cost,'time','cost [$ W^-1]',6)
        % catchability
        plot_domain1D(1:length(catchability),catchability,'time','catchability [m^2 W^-1 s^-1]',7)
    end
    % Save forcing
    Economical.price=price;
    Economical.cost=cost;
    Economical.catchability=catchability;
    save('Economical.mat','Economical','-v7.3')
end

% Load and convert regulation forcing ************  
if create_regulation
    disp('Create regulation forcing')
    % Create effort target ***********************
    efftarg=get_var(efftarg_path,efftarg_var,efftarg_dim);
    % Create social enforcement ******************
    % GOOD LUCK ! (MODIF)
    % Plot forcings ******************************
    if plot_input
        % Efftarg
        plot_domain2D(squeeze(efftarg(:,:,1,1)),lon,lat,mask,'Effort target [??]',8)
    end
    % Save forcing
    Regulation.efftarg=efftarg;
    %(MODIF)
    save('Regulation.mat','Regulation','-v7.3')
end
disp('Done !')



