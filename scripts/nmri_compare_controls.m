function [ all_subjects ] = nmri_statistics( all_subjects, opt, params )
%[ all_subjects ] = nmri_statistics( all_subjects, opt, params )
%   Will do generate a tabulated overview of metrics (all that are readable)
%   via nmri_load_metrics for any number of subject
%
% all_subjects      = cell array of subjects (see nmri_all_subject)
% opt               = struct array of options
%  .report_dir      = target directory for reports
%  .report_regional = include regional information (SUMA ROIs), true/false (logical)
%  .metrics         = cell array of metrics to use
%                     will take from params (+ power, power_noise), if not set
%  .freqs           = cell array of frequencies to use (either char name or index)
%                     will take all from params, if not set
%  .abs             = take the absolute of the metric or not, true/false (logical)
%                     default, true
%  .compare_cx      = compare against all controls (excluding this one, in case)
%  .compare_cx_filt = filter to find controls in all_subjects (see nmri_filter_subjects)
%           default : .id='C\..*'

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
 opt.report=true;
 opt.report_regional=true;
end

if ~isfield(opt,'report') || ~islogical(opt.report)
 opt.report=false;
else
 % we want reports
 if ~isfield(opt,'report')
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

if ~isfield(opt,'abs') || ~islogical(opt.abs)
 opt.abs=true;
end
if opt.abs
 abs_txt='abs';
else
 abs_txt='non_abs';
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

if ~isfield(opt,'compare_cx') || ~islogical(opt.compare_cx)
 opt.compare_cx=false;
else
 if ~isfield(opt,'compare_cx_filt') || ~ischar(opt.compare_cx_filt)
  opt.compare_cx_filt.id='C\..*';
 end
end

%% Setup vars
if opt.report
 stat_cols=6; % always id, exam, metric, freq, gobal_avg
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
 stat_out(1,1:6)={'Subject_ID','Exam_ID','Metric','Freq','Abs','Gobal_Avg'};
 if opt.report_regional
  stat_out(1,7:end)=use_rois_txt';
 end
 stat_row=2;
end


% Main Loop - Per Metric and Per Frequency
for m=1:length(opt.metrics)
 for f=1:length(opt.metrics)
  
  % get normative range, if needed
  if opt.compare_cx
   [cx_metrics]=nmri_load_metrics(nmri_filter_subjects(all_subjects,opt.compare_cx_filt),opt.metrics{m},freq_n(f),params);
  end

  % loop over subjects
  for n=1:length(all_subjects)
   
   % Load the data - for the current subject
   [metric, suma, subjects_back]=nmri_load_metrics(all_subjects(n),opt.metrics{m},freq_n(f),params);
   if ~isequal(subjects_back,all_subjects(n))
    warning(['Could not find the metric ' opt.metrics{m} '/' freq_txt{f} ' for ' all_subjects{n}.id ' -- skipping this one'])
    continue
   end
   
   if isempty(ref_suma)
    ref_suma=suma;
   else
    % check for SUMA annotation equality
    if opt.report_regional && ~isequal(ref_suma{1}.annot_key,suma{1}.annot_key);
     error(['There is a mismatch of SUMA annotations for ' all_subjects{n}.id ' and the first subject - cannot do regional analyis'])
    end
   end
    
   % deal with abs
   if opt.abs
    metric{1}=abs(metric{1});
   end
   
   % Reports
   if opt.report
    % make average overall
    all_mean=metric{1};
    while numel(all_mean)>1
     all_mean=nanmean(all_mean); 
    end
    % put data
    stat_out{stat_row,1}=all_subjects{n}.id;
    stat_out{stat_row,2}=all_subjects{n}.exam_id;
    stat_out{stat_row,3}=opt.metrics{m};
    stat_out{stat_row,4}=freq_txt{f}; 
    stat_out{stat_row,5}=abs_txt; 
    stat_out{stat_row,6}=all_mean; 
    if opt.report_regional
     % now parse the regions
     for r=1:length(use_rois_idx)
      all_mean=metric{1}(suma{1}.annot==use_rois_idx(r));
      while numel(all_mean)>1
       all_mean=nanmean(all_mean); 
      end
      stat_out{stat_row,6+r}=all_mean; 
     end
    end
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
 % Add Java POI Libs to matlab javapath
 javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-3.8-20120326.jar'));
 javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-ooxml-3.8-20120326.jar'));
 javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/poi-ooxml-schemas-3.8-20120326.jar'));
 javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/xmlbeans-2.3.0.jar'));
 javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/dom4j-1.6.1.jar'));
 javaaddpath(fullfile(analysis_dir,'scripts','utilities','xlwrite','poi_library/stax-api-1.0.1.jar'));
 % Generate XLSX file
 xlwrite(fullfile(opt.report_dir,['metric_summary_' date '.xlsx']), stat_out, 'metric_summary');
end



end

