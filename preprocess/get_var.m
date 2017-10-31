%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Load datasets and variables and align dimensions to dim
%**************************************************************************
function var = get_var(path,var,dim)

    % Exit if user defined forcing
    if strcmp(var,'udef')
        return
    end

    % Load dataset
    if exist(path,'file')
        load(path);
    else
        disp('Hum, double-check this path:')
        disp(path)
    end
    
    % Extract variable
    varname=var;
    [struct field]=strtok(var,'.');
    if exist(struct,'var')
        if isempty(field)
            var=eval(var);
        elseif isfield(eval(struct),field(2:end))
            var=eval(var);
        else
            disp('Aye, double-check this variable:')
            disp(varname)
        end
    else
        disp('Aye, double-check this variable:')
        disp(varname)
    end
    
    % Align and average dimensions
    dimmat=ones(1,5);
    dimmat(1:size(dim,2))=dim;
    dimvar=size(var);
    % lat/lon
    pos1=find(dimvar==dimmat(1));
    pos2=find(dimvar==dimmat(2));
    if isempty(pos1) && isempty(pos2)
        disp('Aye, maybe problem with dimension lat/lon for:') 
        disp(varname)
        return
    else
        if isempty(pos1)
            if pos2==1
                var=repmat(var',[dim(1),1]);
            else
                var=repmat(var,[dim(1),1]);
            end
        end
        if isempty(pos2)
            if pos1==1
                var=repmat(var,[1,dim(2)]);
            else
                var=repmat(var',[1,dim(2)]);
            end
        end
    end
    dimvar=size(var);
    pos1=find(dimvar==dimmat(1));
    pos2=find(dimvar==dimmat(2));
    % group
    pos4=find(dimvar==dimmat(4));
    % ensemble
    pos5=find(dimvar==dimmat(5));
    % Add dimensions here if you wish...
    %pos6=find(dimvar==dimmat(6));
    % time
    time=dimvar(~ismember(dimvar,squeeze([dimmat(1) dimmat(2) dimmat(4) dimmat(5)])));
    if ~isempty(time)
        if mod(time,dimmat(3))==0
            if time==dimmat(3)
                pos3=find(dimvar==dimmat(3));        
            else
                pos3=find(dimvar==time);
                if pos3==3
                    tmpvar=var(:,:,1:dimmat(3))/time*dimmat(3);
                    for t=1:time/dimmat(3)-1;
                        tmpvar=tmpvar+var(:,:,1+t*dimmat(3):dimmat(3)+t*dimmat(3))/time*dimmat(3);
                    end
                else
                    disp('Aye, no case exist for time in this dimension... Implement it for:') 
                    disp(varname)
                end
                var=tmpvar;
            end
        else     
            disp('Eish, check the time dimension for:') 
            disp(varname)
        end
    else
        pos3=[];
    end
    % permute
    var=permute(var,[pos1,pos2,pos3,pos4,pos5]);
    
end