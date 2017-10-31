%**************************************************************************
% BOATS PREPROCESS SUBFUNCTION
% Plot ecological/economical 1D domains  
%**************************************************************************
function plot_domain1D(X,Y,Xname,Yname,numfig)
        
        % Plot parameters
        font=16;
        line=3;
        
        % plot domain
        figure(numfig)
        plot(X,Y)
        xlabel(sprintf('%s',Xname),'fontsize',font)
        ylabel(sprintf('%s',Yname),'fontsize',font)
        set(gca,'FontSize',font,'linewidth',line) 
        box on
        
end
