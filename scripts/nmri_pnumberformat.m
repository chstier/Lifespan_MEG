function [ txt ] = nmri_pnumberformat( num, crit1, crit2 )
%[ txt ] = nmri_numberformat( num )
%   Automatically format a number in p-value logic

if ~exist('crit1','var')
 crit1=0.05;
end

if ~exist('crit2','var')
 crit2=0.5;
end

if num<0.001
 txt='<0.001';
elseif num>crit2
 txt=['>' sprintf('%0.1f',crit2)];
else
 txt=sprintf('%0.3f',num);
end
 
if num<crit1
 txt=[txt '*'];
end

end

