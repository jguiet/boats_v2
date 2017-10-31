%**************************************************************************************************************
% FUNCTION modify_parameters.m
% Changes standard parameter values according to user defined or ensembles
% values
%**************************************************************************************************************
function param = modify_parameters(boats,varargin)

%---------------------------------
% Modify path parameters
 param.path = modify(boats.param.path,varargin);

%---------------------------------
% Modify main parameters
 param.main = modify(boats.param.main,varargin);
 
%---------------------------------
% Modify conversion parameters
 param.conversion = modify(boats.param.conversion,varargin);
 
%---------------------------------
% Modify environmental parameters
 param.environment = modify(boats.param.environment,varargin);
 
%---------------------------------
% Modify ecological parameters
 param.ecology = modify(boats.param.ecology,varargin);
 
%---------------------------------
% Modify economical parameters
 param.economy = modify(boats.param.economy,varargin);

%**************************************************************************************************************
% END FUNCTION


%**************************************************************************************************************
% SUBFUNCTIONS
%**************************************************************************************************************
% modify : replace parameters 
function boats = modify(boats,varargin)
% Change parameter values

 Tin = boats;

 if mod(length(varargin{1}),2)~=0
    error('Incorrect number of arguments - must be pairs of property/value');
 end
 nfields = length(varargin{1})/2;
 varPairs = reshape(varargin{1},2,nfields)';

 Tout = Tin;

 for i=1:nfields
     thisPair = varPairs(i,:);
     thisName = thisPair{1};
     thisVal  = thisPair{2};
     % check if property is a string
     if ~strcmp(class(thisName),'char')
        error(['Property ' num2str(i) ' should be a string']);
     end
     % check field is in Tin
     if isfield(Tin,thisName)
        thisValOrig = Tin.(thisName);
        % check class of field is correct
        if class(thisVal)~=class(thisValOrig)
           error(['In property ' num2str(i) ' value is '  class(thisVal) ' kind - should be ' class(thisValOrig)]);
        end
        %substitute value
        Tout.(thisName) = thisVal;
     end
  end

 boats = Tout;
