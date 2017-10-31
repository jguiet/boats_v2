%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Arrange the mask with wet nodes = 0 and dry nodes =1  
%**************************************************************************
function mask = arrange_mask(mask,wet,dry)

    masktmp=mask;

    % Values in mask
    val=unique(mask);
    val(isnan(val(1:end-1)))=[];
    
    % Local function
    ismembernan = @(a,b) ismember(a,b) | (isnan(a) & any(isnan(b)));
    
    % Change wet to 0
    if ismembernan(wet,val)
        if isnan(wet)
            masktmp(isnan(mask))=0;
        else
            masktmp(find(mask==wet))=0;
        end
    else
        disp('Problem with wet cells in mask, the value:') 
        disp(wet)
        disp('is not in what i loaded, i only find the values:')
        disp(val)
    end
    
    % Change dry to 1
    if ismembernan(dry,val)
        if isnan(dry)      
            masktmp(isnan(mask))=1;
        else
            masktmp(find(mask==dry))=1;    
        end
    else
        disp('Problem with dry cells in mask, the value:') 
        disp(dry)
        disp('is not in what i loaded, i only find the values:')
        disp(val)
    end
    
    % Update
    mask=masktmp;
    
end