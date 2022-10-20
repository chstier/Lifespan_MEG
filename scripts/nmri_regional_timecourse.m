function [ data , regional_labels] = nmri_regional_timecourse( data, ref_surf, opt)
%[ regional_metrics ] = nmri_regional_values( full_metrics, ref_surf, opt )
%   Will resample a full brain e.g. SUMA surface (1 value/vertex) metric
%   into a regional metric (1 value / region)
%   the regional information will be taken from the ref_surf struct, that
%   could be the all_suma surface
%   
%
% full_metrics      = cell array of metrics, e.g. from nmri_export_metrics
%
% ref_surf          = reference surface struct, e.g. all_suma, needs to
%                     find .annot_key, .annot fields there  
%
% opt               = struct array of options
%  .dim_reduction   = char of statisitcal measures to use for
%                     dimension reduction
%                     default, 'mean' 
%                     available: mean, median, max, min
%
% regional_metrics  = output, regional metric
 

% check opt and set defaults
if ~exist('opt','var')
 % defaults
 opt=[];
end
if ~isfield(opt,'dim_reduction') || ~ischar(opt.dim_reduction) || isempty(opt.dim_reduction)
 opt.dim_reduction='mean';
end

% check surface
if ~exist('ref_surf','var') || ~isstruct(ref_surf) || ~isfield(ref_surf,'annot') || ~isfield(ref_surf,'annot_key') 
 error('The reference surface does not have the required format')
end

% check metric
if ~exist('data','var') || ~isstruct(data)
 error('You need to provide the data as Fieldtrip struct')
end

% check dimensions
Nchan=size(data.trial{1,1},1);
Ntrials=length(data.trial);

Nsuma=length(ref_surf.annot);
if Nchan~=Nsuma
 error(['There is a mismatch of surface vertices and channels'])
end

% determine unique ROIS
use_rois=ref_surf.annot_key{1}>0;
use_rois_txt=ref_surf.annot_key{2}(use_rois);
use_rois_idx=ref_surf.annot_key{1}(use_rois);

regional_trials=cell(1,Ntrials);

%loop per trial
for n=1:Ntrials
 Tlength=size(data.trial{n},2);
 regional_trials{n}=zeros(length(use_rois_idx),Tlength);
 % now parse the regions
 for r=1:length(use_rois_idx)
  % note: lh and rh have the same code in SUMA - we need to consider this
  % differently
  this_roi=ref_surf.annot==use_rois_idx(r);
  if regexp(use_rois_txt{r},'.*_lh_.*')
   this_roi=this_roi & (ref_surf.hemi==1);
  elseif regexp(use_rois_txt{r},'.*_rh_.*')
   this_roi=this_roi & (ref_surf.hemi==2);
  end
 
  if sum(this_roi)>0
   regional_trials{n}(r,1:Tlength)=calc_stat(data.trial{n}(this_roi,:),opt.dim_reduction); 
  end
 end
end
% export region labels
regional_labels=use_rois_txt;

% build new data struct
data.trials=regional_trials;
clear regional_trials
data.label=use_rois_txt;


end


%% Local functions
function val=calc_stat(metric,stat)
 if iscell(metric)
  metric=metric{1};
 end
 switch(stat)
  case 'abs_mean'
   metric=abs(metric);
   val=nanmean(metric,1);
  case 'mean'
   val=nanmean(metric,1);
  case 'pos_mean'
   metric=metric(metric>0);
   val=nanmean(metric,1);
  case 'neg_mean'
   metric=metric(metric<0);
   val=nanmean(metric,1);
  case 'abs_median'
   metric=abs(metric);
   val=median(metric,'omitnan');
  case 'median'
   val=median(metric,'omitnan');
  case 'pos_median'
   metric=metric(metric>0);
   val=median(metric,'omitnan');
  case 'neg_median'
   metric=metric(metric<0);
   val=median(metric,'omitnan');
  case 'min'
   val=nanmin(metric,[],1);
  case 'max'
   val=nanmax(metric,[],1);
  otherwise
   warning(['Unkown statistical operation =' stat])
   val=NaN;
 end
end
