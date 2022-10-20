function [ p, s, levp , c, mtbl, sumtxt ] = nmri_groupcompare( x , l , bf, use_nonparam)
%[ p, s, levp , mtbl] = nmri_groupcompare( x , l , bf)
%   Does (somewhat smart) group comparisons using ANOVA and/or Kruskal
%   Wallis test (n case Levene is not passed)7*2*2*3
%
%
% p = p value
% s = multiple comparisions tabe
% levp = Levene p
% c = multiple comparisons table
% mtbl = multiple comparisions table (corrected and formated)
% bf = additional Bonferroni correction factor (in addition to the multcompare)
% use_nonparam = logical, if true, always use non parametric

if ~exist('use_nonparam','var')
 use_nonparam=0;
end


if ~exist('bf','var')
 bf=1;
end


sumtxt='';

% check equal variance 
levp=vartestn(x,l,'TestType','LeveneAbsolute','Display','off');
fprintf('Levene test for equal variance p=%0.6f\n',levp)
if levp<0.05 || use_nonparam
 fprintf('Equal variance rejected, perform non-parametric test\n')
 [p,~,s]=kruskalwallis(x,l,'off');
 fprintf('Kruskal Wallis test, p=%s\n',nmri_pnumberformat(p))
 fprintf('Rank (%d) = %d\n',[[1:length(s.meanranks)];s.meanranks])
 means=s.meanranks;
else
 fprintf('Equal variance assumed, perform ANOVA\n')
 [p,~,s]=anova1(x,l,'off');
 fprintf('ANOVA p=%s\n',nmri_pnumberformat(p))
 fprintf('Mean (%d) = %d\n',[[1:length(s.means)];s.means])
 means=s.means;
end

% now do multiple comparision run
[c,~,~]=multcompare(s,'CType','bonferroni','Display','off');



% make a p value matrix
mtbl=cell(size(s.gnames,1)+1);

for i=1:size(c,1)
 %loop over all comparisions
 if isnumeric(s.gnames{c(i,1)})
  txtA=sprintf('%d',s.gnames{c(i,1)});
 else
  txtA=s.gnames{c(i,1)};
 end
 if isnumeric(s.gnames{c(i,2)})
  txtB=sprintf('%d',s.gnames{c(i,2)});
 else
  txtB=s.gnames{c(i,2)};
 end
 mp=nmri_pnumberformat(c(i,6)*bf);
 if c(i,6)*bf<0.05
  if means(c(i,1))<means(c(i,2))
   txtC='<';
  else
   txtC='>';
  end
 else
  txtC='~';
 end
 fprintf('%s %s %s:  p = %s\n',txtA,txtC,txtB,mp)
 sumtxt=[sumtxt sprintf('%s %s %s\t%s\n',txtA,txtC,txtB,mp)];
 
 mtbl{c(i,1)+1,c(i,2)+1}=sprintf('%s (%s%s%s)',mp,txtA,txtC,txtB);
 if isempty(mtbl{c(i,1)+1,1})
  mtbl{c(i,1)+1,1}=txtA;
  mtbl{1,c(i,1)+1}=txtA;
 end
 if isempty(mtbl{1,c(i,2)+1})
  mtbl{1,c(i,2)+1}=txtB;
  mtbl{c(i,2)+1,1}=txtB;
 end
end

% remove unsensible rows/cols
mtbl(end,:)=[]; % last row
mtbl(:,2)=[]; % second col





end

