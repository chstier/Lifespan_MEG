function [pval, realrho, realpval, nullrho] = cs_pvalvsNull_spearman(data1,data2,permno,wsname)
%% This script was build based on the original script by Aaron Alexander-Bloch & Siyuan Liu, pvalvsNull.m, 2018-04-22
% written and edited by Christina Stier, 2022

% Calculate the p-value of correlation between two surface maps based on
% the null distribution of spins of map 1

%load the saved workspace from SpinPermu
load(wsname)

%calculate the real Spearman's correlation between two interested maps
[realrho, realpval] =corr(data1,data2, 'rows','complete', 'type', 'Spearman'); % 'rows','complete' to exclude NaN's
 
% test the observed rho against null described by SpinPermu
nullrho=[];
for i=1:permno
tempdata=cat(2,bigrotl(i,:),bigrotr(i,:))';
nullrho=cat(1,nullrho,corr(tempdata,data2, 'rows','complete', 'type', 'Spearman')); % 'rows','complete' to exclude NaN's
end
%assuming sign is preserved, calculate the probability that the observed
%correlation coeffcient is above the null distribution
pval=length(find(abs(nullrho)>abs(realrho)))/permno; % added abs() 07/31/2020



