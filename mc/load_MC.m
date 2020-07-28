%**************************************************************************************************************
% FUNCTION set_MC.m
% Determine Monte Carlo random parameters for the simulation, from rangesc
%**************************************************************************************************************
function ens_param = load_MC(MParams,Community,replicateid,ens_param)

% Distinct indexes for Pelagic and Demersal communities
if strcmp(Community,'Pelagic')
    id = 1;
elseif strcmp(Community,'Demersal')
    id = 2;
end

%--------------------------------------------------------
% Update variable depending on range and distribution
for indp=1:size(MParams,1)
   ParName = MParams{indp,1};
   ProbFun = lower(MParams{indp,7});
   ParMean = MParams{indp,2};
   ParStd  = MParams{indp,3};
   %--------------------------------------------------------
   rndgen = rand(1);
   switch lower(ProbFun)
   case {'uniform'}
       ParLow = ParMean - sqrt(3)*ParStd;
       ParHig = ParMean + sqrt(3)*ParStd;
       ens_param.(genvarname(ParName))(replicateid,id) = ParLow + rndgen * (ParHig-ParLow);
   case {'gauss','gaussian','norm','normal'}
       ens_param.(genvarname(ParName))(replicateid,id) = norminv(rndgen,ParMean,ParStd);
   case {'gamma'}
       ParShape = ParMean.^2/ParStd.^2;
       ParScale = ParStd.^2/ParMean;
       ens_param.(genvarname(ParName))(replicateid,id) = gaminv(rndgen,ParShape,ParScale);
   case {'lognormal'}
       ens_param.(genvarname(ParName))(replicateid,id) = logninv(rndgen,ParMean,ParStd);
   case {'ens_param.tau0'}
       ParLow = (ens_param.tau0-ParStd)/(log10(MParams{indp,4}/MParams{indp,6}));
       ParHig = (ens_param.tau0-ParMean)/(log10(MParams{indp,5}/MParams{indp,6}));
       ens_param.(genvarname(ParName))(replicateid,id) = - rndgen * (min(ParLow,ParHig));
   otherwise
       error([ProbFun ' distribution not implemented']);
   end
   % Duplicate Pelagic Community for the Demersal community in case the
   % parameters are identical
   if strcmp(Community,'Pelagic')
       ens_param.(genvarname(ParName))(replicateid,id+1) = ens_param.(genvarname(ParName))(replicateid,id);
   end
end   % indp - parameter assignement

ThisClock = clock;
ThisDate  = datestr(now,'mmm-dd-yy');
RndNum   = floor(rand(1,1)*(10000));
ens_param.name(replicateid,:) = [ '_' ThisDate '-' num2str(ThisClock(4)) 'h' num2str(ThisClock(5)) 'm' '-r' num2str(RndNum,'%04d')];

%**************************************************************************************************************
% END FUNCTION

