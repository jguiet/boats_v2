function [corr_P_H2F pval_P_H2F corr_S_H2F pval_S_H2F] = nancorr(x,y)

 validdata1 = ((x ~= 0) & (~isnan(x)));
 validdata2 = ((y ~= 0) & (~isnan(y)));
 validdataAll = validdata1 & validdata2;

 [corr_P_H2F pval_P_H2F] = corr(x(validdataAll),y(validdataAll),'type','Pearson');
 [corr_S_H2F pval_S_H2F] = corr(x(validdataAll),y(validdataAll),'type','Spearman');
 
end
