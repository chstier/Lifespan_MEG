function [ all_subjects ] = nmri_compare_patient( all_subjects, opt, params )
%[ subject ] = nmri_compare_patient( subject )
%   will compare a single patient against the available controls
%   will automatically match by exam_id and headmodel 

% opt               = struct array of options
%  .report_regional = include regional information (SUMA ROIs), true/false (logical)
%  .metrics         = cell array of metrics to use
%                     will take from params (+ power, power_noise), if not set
%  .freqs           = cell array of frequencies to use (either char name or index)
%                     will take all from params, if not set
%  .dim_reduction   = cell array of statisitcal measures to use for
%                     dimension reduction
%                     default, {'mean'}, i.e. mean 
%                     available: mean, median
%  .global          = cell array of global value modification
%                     'abs': absolute, 'pos': only postive, 'neg': only
%                     negative, '' or 'none':nothing
%  .scale           = cell array of global scaling
%                     'mean', 'median', 'rms': root-mean-square, 
%                     '' or 'none':nothing
%  .alpha           = cell array of alpha/p-value threshold (numeric) and
%                     if to use FDR correction or not (logical). Can be
%                     multiple rows.
%                     default: {0.05,true;0.001,false} = both 0.05 FDR
%                     corr, and 0.001 uncorrected
%  .output          = dir for output
%                     default: <analysis_dir>/<subject_id>/results/<exam_id>
%  .show_raw        = logical, also show/safe raw map (after scaling/etc..)
%                     default: true
%  

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


% check opt and set defaults
if ~exist('opt','var')
 % defaults
 opt=[];
 opt.report_regional=false;
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

if ~isfield(opt,'abs') || ~iscell(opt.abs) || isempty(opt.abs)
 opt.abs={'abs'};
end

if ~isfield(opt,'show_raw') || ~islogical(opt.show_raw) || isempty(opt.show_raw)
 opt.show_raw=true;
end


if ~isfield(opt,'dim_reduction') || ~iscell(opt.dim_reduction) || isempty(opt.dim_reduction)
 opt.dim_reduction={'mean'};
end

if ~isfield(opt,'scale') || ~iscell(opt.scale) || isempty(opt.scale)
 opt.scale={'none','mean'};
end

if ~isfield(opt,'global') || ~iscell(opt.global) || isempty(opt.global)
 opt.global={'abs','none'};
end

if ~isfield(opt,'output') || ~ischar(opt.output) || isempty(opt.output)
 opt.output='<analysis_dir>/<subject_id>/results/<exam_id>';
end

if ~isfield(opt,'alpha') || ~iscell(opt.alpha) || isempty(opt.alpha) || size(opt.alpha,2)~=2
 opt.alpha={0.05,true;0.001,false};
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


if isstruct(all_subjects)
 % in case a single subject is given
 all_subjects={all_subjects};
end



% loop over subjects
for n=1:length(all_subjects)
 % process active subject
 subject=all_subjects{n}; 
 % get modality params
 [ this_params ] = nmri_get_modality_params( params, subject.dtype );
 % do some checks
 if ~isfield(subject,'stats') || ~exist(subject.stats,'file') 
  warning(['Subject ID=' subject.id 'does not seem to be processed - skipping'])
  continue
 end
 if ~strcmp(subject.id(1),'P')
  warning(['Subject ID=' subject.id 'does not seem to be a patient - skipping'])
  continue
 end
 % Get the right SUMA template for this one
 suma_ld=subject.stamps.(['processing_' this_params.headmodel]).params.SUMA_ld;
 suma_fsaverage=fullfile(subject.analysis_dir,'conf',['suma-all-fsaverage-' suma_ld '.mat']);
 if ~exist(suma_fsaverage,'file')
  error(['Could not find the pre-processed SUMA fsaverage for LD= ' suma_ld ', file=' suma_fsaverage])
 end
 load(suma_fsaverage,'suma_all')

 
 % Main Loop - Per Metric and Per Frequency
 for m=1:length(opt.metrics)
  for f=1:length(opt.freqs)
   fprintf('Subject=%s, Exam_ID=%s, Metric=%s\n',subject.id,subject.exam_id,opt.metrics{m})

   % find matching controls
   filters=[];
   filters.exam_id=['^' subject.exam_id '.*'];
   filters.id='C\..*';
   filters.stamps.(['processing_' this_params.headmodel])=true;
   fprintf('Searching for matching controls...\n')
   [ filtered_subjects ] = nmri_filter_subjects ( nmri_all_subjects(subject.analysis_dir), filters );
   fprintf('...%d controls found\n',length(filtered_subjects))
   
   % Load the data - for the current subject and controls
   this_subjects=[{subject}; filtered_subjects];
   [all_metrics, all_suma, subjects_back]=nmri_load_metrics(this_subjects,opt.metrics{m},freq_n(f),params);
   if ~isequal(subjects_back,this_subjects)
    warning(['Could not find the metric ' opt.metrics{m} '/' freq_txt{f} ' for all subjects / controls -- skipping'])
    continue
   end
   % get all suma msk, start with fsaverage
   all_msk=suma_all.msk;
   for si=1:length(all_suma)
    if isequal(size(all_msk),size(all_suma{si}.msk))
     all_msk=all_msk&all_suma{si}.msk;
    else
     error(['Problem with the SUMA vertex count for ID=' filtered_subjects{si}.id])
    end
   end
   % remove from all
   for mi=1:length(all_metrics)
    if size(all_metrics{mi},2)>1 && size(all_metrics{mi},1)>1
     % two dims
     all_metrics{mi}(~all_msk,~all_msk)=NaN;  
    else
     % one dim
     all_metrics{mi}(~all_msk)=NaN;   
    end
   end
   
   % now make fieldtrip compatible subject structs based on fsaverage
   % first take from the template source model
   sourceStruct            = [];
   sourceStruct.pos        = suma_all.pos;
   sourceStruct.inside     = all_msk;
   sourceStruct.outside    = ~all_msk;

   allStructs = cell(length(this_subjects),1);
   % loop through subjects and create a source structure for each
   for iSubject = 1:length(this_subjects)
    % store it as a structure
    allStructs{iSubject} = sourceStruct;    
   end
  
   
   % loop over all dim_reductions
   for dim_i=1:length(opt.dim_reduction)
    % now check and dim reduce each metric
    dim_metrics=all_metrics;
    for mi=1:length(dim_metrics)
     dim_metrics{mi}=reduce_dim(dim_metrics{mi},opt.dim_reduction{dim_i});
     % and rotate if needed
     if (size(dim_metrics{mi},2)>size(dim_metrics{mi},1))
      dim_metrics{mi}=dim_metrics{mi}';
     end
    end 
    
    % now loop over global abs/pos/neg
    for global_i=1:length(opt.global)
     % now apply global 
     global_metrics=dim_metrics;
     for mi=1:length(global_metrics)
      switch(opt.global{global_i})
      case 'abs'
       global_metrics{mi}=abs(global_metrics{mi});
      case 'pos'
       global_metrics{mi}(global_metrics{mi}<0)=NaN;
      case 'neg'
       global_metrics{mi}(global_metrics{mi}>0)=NaN;
      end
     end
     
     % now loop over scale
     for scale_i=1:length(opt.scale)
      scale_metrics=global_metrics;
      for mi=1:length(scale_metrics)
       switch(opt.scale{scale_i})
        case 'mean'
         val=nanmean(reshape(scale_metrics{mi},numel(scale_metrics{mi}),1));
        case 'median'
         val=median(reshape(scale_metrics{mi},numel(scale_metrics{mi}),1),'omitnan');
         case 'rms'
          val=rms(reshape(scale_metrics{mi},numel(scale_metrics{mi}),1));
        otherwise
         val=1;
       end
       scale_metrics{mi}=scale_metrics{mi}/val;
      end
      
      %% Final stats stage now
      % Note that if you want to run a two-sided test,
      % you have to split the critical alpha value by
      % setting cfg.correcttail = 'alpha'; i.e. this sets cfg.alpha = 0.025,
      % corresponding to a false alarm rate of 0.05 in a two-sided test
      
      for alpha_i=1:size(opt.alpha,1)
       cfg                  = [];
       cfg.method           = 'analytic';
       cfg.statistic        = 'nmri_statfun_indepsamplesT';
       cfg.alpha            = opt.alpha{alpha_i,1};
       %cfg.correcttail      = 'alpha';
       cfg.parameter        = 'param';
       if opt.alpha{alpha_i,2} % check for FDR
        cfg.correctm         = 'fdr';
       end

       % make the design - very special format in fieldtrip, working with
       % rows
       cfg.design=zeros(1,length(this_subjects)); 
       for des_i=1:length(this_subjects)
        if strcmp(this_subjects{des_i}.id(1),'C')
         cfg.design(1,des_i)=2;
        elseif strcmp(this_subjects{des_i}.id(1),'P')
         cfg.design(1,des_i)=1;
        else
         error('Not patient or control, should not happen')
        end
       end
       cfg.ivar             = 1; % the 1st row in cfg.design contains the independent variable

       % now feed in the scaled metric
       for iSubject = 1:length(this_subjects)
        allStructs{iSubject}.param =  scale_metrics{iSubject};
       end

       % TBD: smoothing not yet implemented

       stats                = ft_sourcestatistics(cfg,allStructs{:});

       % cont - here, save stats. 
       % save stats and make image
       outputdir=opt.output;
       outputdir=strrep(outputdir,'<subject_id>',subject.id);
       outputdir=strrep(outputdir,'<exam_id>',subject.exam_id);
       outputdir=strrep(outputdir,'<analysis_dir>',subject.analysis_dir);
       if strcmp(opt.global{global_i},'pos') 
        global_txt='pos_only';
       elseif strcmp(opt.global{global_i},'neg')
        global_txt='neg_only';
       elseif strcmp(opt.global{global_i},'abs')
        global_txt='abs';
       else
        global_txt='no_abs';
       end
       if strcmp(opt.scale{scale_i},'mean') 
        scale_txt='mean_scaled';
       elseif strcmp(opt.scale{scale_i},'median')
        scale_txt='median_scaled';
       elseif strcmp(opt.scale{scale_i},'rms')
        scale_txt='rms_scaled';
       else
        scale_txt='not_scaled';
       end
       if strcmp(opt.dim_reduction{dim_i},'mean') 
        dim_txt='';
       else
        dim_txt=['_dimred_' opt.dim_reduction{dim_i}];
       end
       
       outputdir=fullfile(outputdir,[global_txt '_' scale_txt dim_txt]);
       corr_txt=['p' num2str(opt.alpha{alpha_i,1})];
       if opt.alpha{alpha_i,2}
        corr_txt=[corr_txt 'FDR'];
       else
        corr_txt=[corr_txt 'uncorr'];
       end
       if ~exist(outputdir,'dir')
        mkdir(outputdir)
       end
       
       filebase=fullfile(outputdir,[subject.id '_' subject.exam_id '_vs_' num2str(length(this_subjects)-1) 'Cx_' opt.metrics{m} '_' freq_txt{f} '_' corr_txt]);
       save([filebase '.mat'],'stats','this_subjects')
       pcfg=[];
       pcfg.per_hemi=1;
       pcfg.per_cortex=0;
       pcfg.clim=[2 8];
       % Increase 
       pcfg.title=[subject.id '_' subject.exam_id ' vs ' num2str(length(this_subjects)-1) 'Cx - ' opt.metrics{m} ' - ' freq_txt{f} ' - ' corr_txt ' - Increase\n' global_txt ' - ' scale_txt dim_txt];
       func=stats.stat;
       func(~stats.mask)=0;
       func(func<stats.critval(2))=0;
       act_vert=sum(func>0);
       if act_vert>0
        pcfg.title=[pcfg.title '\nSuprathreshold vertices = ' num2str(act_vert)];
       else
        pcfg.title=[pcfg.title '\nNo suprathreshold vertices found'];
       end
       pcfg.opathresh=stats.critval(2)*2;
       pcfg.colormap='autumn';
       pcfg.output=[filebase '_increase.png'];
       hFig=nmri_plot_surface_suma(all_suma{1},func,pcfg);
       close(hFig)
       % Decrease 
       pcfg.title=[subject.id '_' subject.exam_id ' vs ' num2str(length(this_subjects)-1) 'Cx - ' opt.metrics{m} ' - ' freq_txt{f} ' - ' corr_txt ' - Decrease\n' global_txt ' - ' scale_txt dim_txt];
       func=-stats.stat;
       func(~stats.mask)=0;
       func(func<-stats.critval(1))=0;
       act_vert=sum(func>0);
       if act_vert>0
        pcfg.title=[pcfg.title '\nSuprathreshold vertices = ' num2str(act_vert)];
       else
        pcfg.title=[pcfg.title '\nNo suprathreshold vertices found'];
       end
       pcfg.opathreshthresh=-(stats.critval(1)*2);
       pcfg.colormap='cool';
       pcfg.output=[filebase '_decrease.png'];
       hFig=nmri_plot_surface_suma(all_suma{1},func,pcfg);
       close(hFig)
       
      end % Alpha loop
      
      % print raw
      if opt.show_raw
       filebase=fullfile(outputdir,[subject.id '_' subject.exam_id '_' opt.metrics{m} '_' freq_txt{f} '_raw']);
       save([filebase '.mat'],'scale_metrics','this_subjects')

       pcfg=[];
       pcfg.per_hemi=1;
       pcfg.per_cortex=0;
       pcfg.title=[subject.id '_' subject.exam_id ' - Raw - ' opt.metrics{m} ' - ' freq_txt{f} '\n' global_txt ' - ' scale_txt dim_txt];
       func=scale_metrics{1};
       pcfg.output=[filebase '.png'];
       hFig=nmri_plot_surface_suma(all_suma{1},func,pcfg);
       close(hFig)
      end
      
     end 
    end 
   end 
  end 
 end 
end 


%% Local functions
function metric=reduce_dim(metric,stat)
 if iscell(metric)
  metric=metric{1};
 end
 if size(metric,2)>1 && size(metric,1)>1
  % reduce
  switch(stat)
   case 'abs_mean'
    metric=abs(metric);
    metric=nanmean(metric,1);
   case 'mean'
    metric=nanmean(metric,1);
   case 'pos_mean'
    metric=metric(metric>0);
    metric=nanmean(metric,1);
   case 'neg_mean'
    metric=metric(metric<0);
    metric=nanmean(metric,1);
   case 'abs_median'
    metric=abs(metric);
    metric=median(metric,'omitnan');
   case 'median'
    metric=median(metric,'omitnan');
   case 'pos_median'
    metric=metric(metric>0);
    metric=median(metric,'omitnan');
   case 'neg_median'
    metric=metric(metric<0);
    metric=median(metric,'omitnan');
   otherwise
    warning(['Unkown statistical operation =' stat])
    metric=NaN(1,size(metric,1));
  end
 end
end


end