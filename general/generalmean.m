%-----------------------------------------------------------------------------------------
% function generalmean(varin,tdim,istart,iend)
% Calculate temporal mean of a variable
% Used by the boats0d_time_average.m function
%-----------------------------------------------------------------------------------------

function varout = generalmean(varin,tdim,istart,iend)

    tsize = size(varin);
    ndim = length(tsize);

    if ndim < tdim
       error(['Specified dimension is larger than array dimension']);
    end 

    if tsize(tdim)<istart | tsize(tdim)<iend
       error(['Indexes for averaging are out of bound']);
    end 

    order = [tdim setdiff(1:ndim,tdim)];
    varin = permute(varin, order);
    tsize = size(varin);
    varin = reshape(varin, tsize(1), []);
    varin = varin(istart:iend,:);
    varin = reshape(varin, [length(istart:iend) tsize(2:end)]);
    varout = nanmean(varin, 1);
    varout = permute(varout, [2:ndim 1]);

end

%-----------------------------------------------------------------------------------------
% END OF SCRIPT