%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Create basic user defined forcing
%**************************************************************************
function var = udef_var(type, dim, ref1, ref2, ref3, ref4)

    % Exit if no user defined type
    if strcmp(type,'')
        return
    end

    % Constant forcing
    if strcmp(type,'cst')
        var=squeeze(ones(dim))*ref1;
    end
    
    % Constant rate increase
    if strcmp(type,'rate')
        for t=1:dim(3)
            if t==1
                var(1:dim(1),1:dim(2),t)=ones(dim(1),dim(2))*ref1*(1+ref2)^(1/12);
            else
                var(1:dim(1),1:dim(2),t)=var(1:dim(1),1:dim(2),t-1)*(1+ref2)^(1/12);
            end
        end
        var=squeeze(var);
    end
    
    
    % Undefined forcing type
    if strcmp(type,'??')
        disp('Aye, undefined forcing type, ask Jerome')
        return
    end
end