function [corr_P pval_P corr_S pval_S] = nancorrlog10(x,y)

 validdata1 = ((x ~= 0) & (~isnan(x)));
 validdata2 = ((y ~= 0) & (~isnan(y)));
 validdata3 = (x./y < 10000);
 validdataAll = validdata1 & validdata2 & validdata3;

 if (nansum(validdataAll) > 5)
 
   xA = log10(x(validdataAll));
   yA = log10(y(validdataAll));

   [corr_P pval_P] = corr(xA,yA,'type','Pearson');
   [corr_S pval_S] = corr(xA,yA,'type','Spearman');
 
 else

   corr_P = NaN;
   pval_P = NaN;
   corr_S = NaN;
   pval_S = NaN;

 end

end
