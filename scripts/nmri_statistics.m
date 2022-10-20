function [ stat_out ] = nmri_statistics( all_subjects, opt, params )
%[ all_subjects ] = nmri_statistics( all_subjects, opt, params )
%   Will do generate a tabulated overview of metrics (all that are readable)
%   via nmri_load_metrics for any number of subject
%
% all_subjects      = cell array of subjects (see nmri_all_subject)
%
% opt               = struct array of options
%  .report_dir      = target directory for reports
%  .report_regional = include regional information (SUMA ROIs), true/false (logical)
%  .metrics         = cell array of metrics to use
%                     will take from params (+ power, power_noise), if not set
%  .freqs           = cell array of frequencies to use (either char name or index)
%                     will take all from params, if not set
%  .statistics      = cell array of statisitcal measures to use
%                     default, {'abs_mean'}, i.e. mean of absolute values
%                     available: abs_mean, mean, abs_median, median (with our without abs) 
%                     pos_mean, neg_mean, pos_median, neg_median (only >0 or <0)
%
% stat_out          = cell array table with values and captions


% check the params struct
if (~exist('params','var'))
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in the current path, or have it in the call ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the paramter file (analysis_params.m)')  
  end
 end
end


if isstruct(all_subjects)
 % in case a single subject is given
 all_subjects{1}=all_subjects;
end


% check opt and set defaults
if ~exist('opt','var')
 % defaults
 opt.report_regional=true;
end

% we always report
opt.report=true;

if ~isfield(opt,'report') || ~islogical(opt.report)
 opt.report=false;
else
 % we want reports
 if ~isfield(opt,'report_dir')
  opt.report_dir=fullfile(pwd,'reports');
 end
 if ~isfield(opt,'report_regional') || ~islogical(opt.report_regional)
  opt.report_regional=false;
 end
end

if ~isfield(opt,'metrics') || ~iscell(opt.metrics) || isempty(opt.metrics)
 % as a default we take everything
 opt.metrics={};
 all_metrics=[params.con_method {'power','power_noise'}];
 for i=1:length(all_metrics)
  if strcmp(all_metrics{i},'coh')
   opt.metrics{end+1}=[all_metrics{i} '_real'];
   opt.metrics{end+1}=[all_metrics{i} '_img'];
  else
   opt.metrics{end+1}=all_metrics{i};
  end
 end
end

if ~isfield(opt,'statistics') || ~iscell(opt.statistics) || isempty(opt.statistics)
 opt.statistics={'abs_mean','mean'};
end


if ~isfield(opt,'freqs') || ~iscell(opt.freqs) || isempty(opt.freqs)
 % as a default we take everything
 opt.freqs=1:length(params.freqs);
end
% check freq text or number
freq_txt=cell(length(opt.freqs),1);
freq_n=zeros(length(opt.freqs),1);
for i=1:length(opt.freqs)
 if isnumeric(opt.freqs(i))
  freq_txt(i)=params.freqsNames(opt.freqs(i));
  freq_n(i)=opt.freqs(i);
 elseif iscell(opt.freqs(i)) && ischar(opt.freqs{i})
  freq_txt{i}=opt.freqs{i};
  freq_n=find(strcmp(params.freqsNames,opt.freqs{i}));
 else
  error('Could not determine frequency requested')
 end
end
%% Setup vars

stat_items={'Subject_ID','Exam_ID','Metric','Freq','Stat','Global'};
stat_cols=length(stat_items); 
 
% load SUMA as ref
[~, ref_suma, ~]=nmri_load_metrics(all_subjects(1),opt.metrics{1},freq_n(1),params);
if opt.report_regional
 % determine unique ROIS
 use_rois=ref_suma{1}.annot_key{1}>0;
 use_rois_txt=ref_suma{1}.annot_key{2}(use_rois);
 use_rois_idx=ref_suma{1}.annot_key{1}(use_rois);
 stat_cols=stat_cols+sum(use_rois);
end
stat_out=cell((length(all_subjects)*length(opt.metrics)*length(opt.metrics))+1,stat_cols);

% Print Header
stat_out(1,1:length(stat_items))=stat_items;
if opt.report_regional
 stat_out(1,7:end)=use_rois_txt';
end
stat_row=2;


% Main Loop - Per Metric and Per Frequency
for m=1:length(opt.metrics)
 for f=1:length(opt.freqs)
  % loop over subjects
  for n=1:length(all_subjects)
   % Load the data - for the current subject
   [metric, suma, subjects_back]=nmri_load_metrics(all_subjects(n),opt.metrics{m},freq_n(f),params);
   if ~isequal(subjects_back,all_subjects(n))
    warning(['Could not find the metric ' opt.metrics{m} '/' freq_txt{f} ' for ' all_subjects{n}.id ' -- skipping this one'])
    continue
   end
   
   % check for SUMA annotation equality
   if opt.report_regional && ~isequal(ref_suma{1}.annot_key,suma{1}.annot_key);
    error(['There is a mismatch of SUMA annotations for ' all_subjects{n}.id ' and the first subject - cannot do regional analyis'])
   end
   
   % now loop per stat requested
   for s=1:length(opt.statistics)
      
    % put default stat_items
    for si=1:length(stat_items)
     switch stat_items{si}
      case 'Subject_ID'
       stat_out{stat_row,si}=all_subjects{n}.id;
       case 'Exam_ID'
       stat_out{stat_row,si}=all_subjects{n}.exam_id;
      case 'Metric'
       stat_out{stat_row,si}=opt.metrics{m};
      case 'Freq'
       stat_out{stat_row,si}=freq_txt{f};       
      case 'Stat'
       stat_out{stat_row,si}=opt.statistics{s};
      case 'Global'
       stat_out{stat_row,si}=calc_stat(metric{1},opt.statistics{s});
     end
    end
    
    if opt.report_regional
     % now parse the regions
     for r=1:length(use_rois_idx)
      % note: lh and rh have the same code in SUMA - we need to consider this
      % differently
      this_roi=suma{1}.annot==use_rois_idx(r);
      if regexp(use_rois_txt{r},'.*_lh_.*')
       this_roi=this_roi & (suma{1}.hemi==1);
      elseif regexp(use_rois_txt{r},'.*_rh_.*')
       this_roi=this_roi & (suma{1}.hemi==2);
      end
      % check if 1D or 2D
      if size(metric{1},1)==1 || size(metric{1},2)==1
       %1D
       stat_out{stat_row,si+r}=calc_stat(metric{1}(this_roi),opt.statistics{s}); 
      else
       %2D
       stat_out{stat_row,si+r}=calc_stat(metric{1}(this_roi,this_roi),opt.statistics{s}); 
      end
     end
    end
    
    % next subject
    stat_row=stat_row+1;
   end
  end
 end
end


% save report
if opt.report
 if ~exist(opt.report_dir,'dir')
  mkdir(opt.report_dir)
 end
  
 % save as csv
 nf_csvwrite(fullfile(opt.report_dir,['metric_summary_' date '.csv']),stat_out);
 
 % and mat
 save(fullfile(opt.report_dir,['metric_summary_' date '.mat']),'stat_out');
 
 % Initialisation of POI Libs
 % determine analyis dir
 analysis_dir=all_subjects{1}.analysis_dir;
 if ~exist(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-3.8-20120326.jar'),'file')
  analysis_dir=pwd;
 end
 if exist(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-3.8-20120326.jar'),'file')
  % Add Java POI Libs to matlab javapath
  javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-3.8-20120326.jar'));
  javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-ooxml-3.8-20120326.jar'));
  javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-ooxml-schemas-3.8-20120326.jar'));
  javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/xmlbeans-2.3.0.jar'));
  javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/dom4j-1.6.1.jar'));
  javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/stax-api-1.0.1.jar'));
  % Generate XLSX file
  xlwrite(fullfile(opt.report_dir,['metric_summary_' date '.xlsx']), stat_out, 'metric_summary');
 else
  warning('Could not locate xlwrite utility, skipping XLSX generation')
 end
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
  otherwise
   warning(['Unkown statistical operation =' stat])
   val=NaN;
 end
end

end

