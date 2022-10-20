function [ regional_metrics , regional_labels] = nmri_regional_values( full_metrics, ref_surf, opt)
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
if ~exist('full_metrics','var')
 error('You need to provide the full brain metric')
end
if isnumeric(full_metrics)
 % make cell in all cases
 tmp=full_metrics;
 full_metrics={};
 full_metrics{1}=tmp;
elseif ~iscell(full_metrics)
 error('Metrics need to be given as cell array (per subject) or numeric array (single subject)')
end

% check dimensions
N=length(full_metrics);
Nsuma=length(ref_surf.annot);
for n=1:N
 if length(full_metrics{n})~=Nsuma
  error(['There is a mis-match of surface vertices and provided metric for n=' num2str(n)])
 end
end

% determine unique ROIS
use_rois=ref_surf.annot_key{1}>0;
use_rois_txt=ref_surf.annot_key{2}(use_rois);
use_rois_idx=ref_surf.annot_key{1}(use_rois);

regional_metrics=cell(N,1);

% check the dimensions
Ndim=0;
Ni=0;
for n=1:N
 elem=size(full_metrics{n});
 Tdim=sum(elem>1);
 if Ndim==0
  Ndim=Tdim;
 elseif (Ndim~=Tdim)
  error('Mismatch of dimension, all subjects need to have the same number of non-singleton dimensions')
 end
 
 % check elements
 Ti=size(full_metrics{n},1);
 if Ni==0
  Ni=Ti;
 elseif (Ni~=Ti)
  error('Mismatch of matrix elements, all subjects need to have the same number of vertices/nodes')
 end
end


% build the region mapping
region_mapping=zeros(Ni,1);
for r=1:length(use_rois_idx)
 % note: lh and rh have the same code in SUMA - we need to consider this
 % differently
 this_roi=ref_surf.annot==use_rois_idx(r);
 if regexp(use_rois_txt{r},'.*_lh_.*')
  this_roi=this_roi & (ref_surf.hemi==1);
 elseif regexp(use_rois_txt{r},'.*_rh_.*')
  this_roi=this_roi & (ref_surf.hemi==2);
 end
 % build the region mapping table
 region_mapping(this_roi,1)=r;
end

%loop subjects

for n=1:N
 regional_metrics{n}=zeros(length(use_rois_idx),1);
 % now parse the regions
 for r=1:length(use_rois_idx)
  % now map the regions based in Ndim
  if Ndim==1
   regional_metrics{n}(r,1)=calc_stat(full_metrics{n}(region_mapping==r),opt.dim_reduction); 
  elseif Ndim==2
   for rr=1:length(use_rois_idx)
    if r==rr
     regional_metrics{n}(r,rr)=NaN;
    else
     regional_metrics{n}(r,rr)=calc_stat(full_metrics{n}(region_mapping==r,region_mapping==rr),opt.dim_reduction); 
    end
   end
  else
   error(['Ndim = ' num2str(Ndim) ' not available'])
  end
 end
end



% export region labels
regional_labels=use_rois_txt;
end


%% Local functions
function val=calc_stat(metric,stat)
 if iscell(metric)
  metric=metric{1};
 end
 switch(stat)
  case 'abs_mean'
   metric=abs(metric);
   val=nanmean(reshape(metric,numel(metric),1));
  case 'mean'
   val=nanmean(reshape(metric,numel(metric),1));
  case 'pos_mean'
   metric=metric(metric>0);
   val=nanmean(reshape(metric,numel(metric),1));
  case 'neg_mean'
   metric=metric(metric<0);
   val=nanmean(reshape(metric,numel(metric),1));
  case 'abs_median'
   metric=abs(metric);
   val=median(reshape(metric,numel(metric),1),'omitnan');
  case 'median'
   val=median(reshape(metric,numel(metric),1),'omitnan');
  case 'pos_median'
   metric=metric(metric>0);
   val=median(reshape(metric,numel(metric),1),'omitnan');
  case 'neg_median'
   metric=metric(metric<0);
   val=median(reshape(metric,numel(metric),1),'omitnan');
  case 'min'
   val=nanmin(reshape(metric,numel(metric),1));
  case 'max'
   val=nanmax(reshape(metric,numel(metric),1));
  otherwise
   warning(['Unkown statistical operation =' stat])
   val=NaN;
 end
end
