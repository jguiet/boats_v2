%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Plot ecological/economical 2D domains 
%**************************************************************************
function plot_domain2D(domain,lon,lat,mask,titlename,numfig)
        
        % Plot parameters
        font=16;
        line=3;
        
        % plot domain
        figure(numfig)
        if length(size(domain))==2
            % 2D
            pcolor(lon,lat,domain.*(1-mask));
        else
            % Or mean 3D
            pcolor(lon,lat,mean(domain,3).*(1-mask));
        end
        xlabel('lon','fontsize',font)
        ylabel('lat','fontsize',font)
        title(sprintf('%s',titlename),'fontsize',font)
        set(gca,'FontSize',font,'linewidth',line) 
        shading flat
        colorbar  
        
end
