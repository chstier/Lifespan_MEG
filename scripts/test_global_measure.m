function [p,h,stat] = test_global_measure(A,B,group1str,group2str,parameterName,freqs)
% where A and B the input data, which have the size Subjects x Frequencies  
% perform nonparametric Mann-Whitney-Wilcoxon test and plot the figures
% based on work of Adham and Yiwen, see there papers. Adopted a bit to work
% with the pipeline workflow

p = ones(size(freqs)); h = ones(size(freqs));

for iFreq = 1:length(freqs)
    % perform Wilcoxon rank sum test, tests the null hypothesis that data
    % in A and B are from continuous distributions with equal medians
    [p(iFreq),h(iFreq),stat(iFreq)] = ranksum(A(:,iFreq),B(:,iFreq)); 
    % this test is equivalent to a Mann-Whitney U test, see below:
    %mwwtest(A(:,iFreq),B(:,iFreq))
end





% option1: use the standard errors as error bars  
sem1 = std(A); % standard error of the mean1 - switched to just std.dev
sem2 = std(B); % standard error of the mean2 - switched to just std.dev

% option2: use the 25% and 75% percentiles of the data vector as error bars
iqr1 = prctile(A,[25 75]);
iqr2 = prctile(B,[25 75]);

% plot median values with error bars
figure
errorbar((1:length(freqs))-0.075,mean(A),-sem1,sem1,'ko','color',[0.3, 0.3, 0.3],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize',8,'MarkerEdgeColor',[0.3, 0.3, 0.3],'LineWidth',1);
hold on
errorbar((1:length(freqs))+0.075,mean(B),-sem2,sem2,'o','color',[0.7, 0.1, 0.1],'MarkerFaceColor',[0.9, 0.5, 0.5],'MarkerSize',8,'MarkerEdgeColor',[0.7, 0.1, 0.1],'LineWidth',1);

%errorbar((1:length(freqs))-0.075,median(A),iqr1(1,:)-median(A),median(A)-iqr1(2,:),'ko','color',[0.3, 0.3, 0.3],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize',8,'MarkerEdgeColor',[0.3, 0.3, 0.3],'LineWidth',1);
%hold on
%errorbar((1:length(freqs))+0.075,median(B),iqr2(1,:)-median(B),median(B)-iqr2(2,:),'o','color',[0.7, 0.1, 0.1],'MarkerFaceColor',[0.9, 0.5, 0.5],'MarkerSize',8,'MarkerEdgeColor',[0.7, 0.1, 0.1],'LineWidth',1);
% title([ parameterName ' per Frequency Band']);
legend(group1str,group2str,'p<0.05') % put legeneds

% fix the ticks
% set(gca,'XTickLabel',{'','Delta','','Theta','','Alpha','','Beta1','','Beta2','','gamma'})
set(gca,'XTickLabel',[{''} freqs],'FontSize',12,'FontWeight','normal')
% set(gca,'YTickLabel',{'','Delta','Theta','Alpha','Beta1','Beta2','gamma'},'FontSize',14)

set(gca, 'Ticklength', [0 0])
theMax = max([A;B]);

plot(find(h),theMax(find(h)),'r*','MarkerSize',12, 'LineWidth',2); % add a red star when the null hypothesis is rejected
ylabel(parameterName);
if any(find(h))
 % add legend for significance
 legend(group1str,group2str,'p<0.05') % put legeneds
end

end