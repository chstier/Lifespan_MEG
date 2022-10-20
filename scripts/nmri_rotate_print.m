function [ hFig ] = nmri_rotate_print( hFig, this_rot, output, title_txt, clipping)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
% get axes
hAx=findall(hFig,'type','axes');
% set export for printer
fig_count=0;
all_plots='';
plot_list={};
set(hFig,'PaperUnits','inches','PaperPosition',[0 0 4 4]) % set square output
for Ri=1:size(this_rot,1)
 fig_count=fig_count+1;
 view(hAx,this_rot(Ri,:))
 this_plot=[output '_subplot_' sprintf('%02d',fig_count) '.png'];
 all_plots=[all_plots ' ' this_plot];
 print(hFig,this_plot,'-dpng','-r200')
 plot_list{end+1}=this_plot;
end
output_all=[output '.png'];
if exist('clipping','var') && length(clipping)==4
 % apply some  clipping to images
 for n=1:length(plot_list)
  cmd=['convert ' plot_list{n} ' -crop ' sprintf('%dx%d+%d+%d',clipping(1),clipping(2),clipping(3),clipping(4)) ' ' plot_list{n}];
  system(cmd);
 end
end
% merge all plots
cmd=['montage -gravity South -mode concatenate -tile ' num2str(fig_count) 'x -background white' all_plots ' ' output_all];
system(cmd);
 
% and annotate 
cmd=['convert ' output_all ' -font /usr/share/fonts/liberation/LiberationSans-Bold.ttf -pointsize ' num2str(36) ' -fill black -background white -gravity North -splice 0x48 -annotate +0+2 "' title_txt '" ' output_all];
system(cmd);
% delete temporary
system(['rm -f ' all_plots]); 
saveas(hFig,[output '.fig'],'fig');


end

