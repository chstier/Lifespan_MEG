function [ output ] = cs_nf_ft_report_regional( suma_all, metric )
% [ output ] = nf_ft_report_regional( suma_all, metric )
%   Will report reginal statistics for a metric (one value per vertex)
%   output is a cell array, could be saved with nf_csvwrite, nf_xlswrite

% CS Aug 2019
% changes to the orignial script (nf_ft_report_regional): take nanmean over vertices for one regions in order to get a value out of it and not only NaN

% determine unique ROIs
use_rois=suma_all.annot_key{1}>0;
use_rois_txt=suma_all.annot_key{2}(use_rois);
use_rois_idx=suma_all.annot_key{1}(use_rois);

output=cell(length(use_rois_idx)+1,4);
output{1,1}='ROI';
output{1,2}='~0 vertices';
output{1,3}='mean';
output{1,4}='mean (~0)';


% now loop over ROIs
for r=1:length(use_rois_idx)
 % note: lh and rh have the same code in SUMA - we need to consider this
 % differently
 this_roi=suma_all.annot==use_rois_idx(r);
 if regexp(use_rois_txt{r},'.*_lh_.*')
  this_roi=this_roi & (suma_all.hemi==1);
 elseif regexp(use_rois_txt{r},'.*_rh_.*')
  this_roi=this_roi & (suma_all.hemi==2);
 end
 
 if sum(this_roi)>0
  % we have at least one vertex left - so count
  txt=use_rois_txt{r};
  % beautify a bit
  txt=strrep(txt,'wm_','');
  txt=strrep(txt,'lh_','Left ');
  txt=strrep(txt,'rh_','Right ');
  txt=strrep(txt,'Left-','Left ');
  txt=strrep(txt,'Right-','Right ');
  txt=strrep(txt,'_',' ');
  output{r+1,1}=txt; % ROI name
  output{r+1,2}=sum(metric(this_roi)>0); % non zero
  output{r+1,3}=nanmean(metric(this_roi)); % mean
  output{r+1,4}=nanmean(nonzeros(metric(this_roi))); % mean non zero   
 end
end

