function [ subject ] = nmri_processing_sensor(subject, data, params)
%[ subject ] = nmri_processing_sensor(subject, data, params)
%  
% This function will do the processing @ sensor level only
% power @sensor level 
% connectivity @sensor level
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 07/2017 

% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

if ~exist('data','var') || isempty(data) 
 % check if we have cleaned dataset (ICA or not)
 if (isfield(params,'useICA_clean') && params.useICA_clean==1)
  if (isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
   input=subject.cleanICA_dataset;
   useICA=true;
  else
   error('ICA cleaned dataset not found - Run nmri_artifactrejection first')
  end
 else
  if (isfield(subject,'clean_dataset') &&  exist(subject.clean_dataset,'file'))
   input=subject.clean_dataset;
   useICA=false;
  else
   error('Cleaned dataset (w/o ICA) not found - Run nmri_artifactrejection first')
  end  
 end
 % retain subject info, if different
 load(input,'data'); 
 if (~exist('data','var') ) 
  error('Could not load data')
 end
end



%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% make dir / files
if (~isfield(subject,'stats_dir'))
 subject.stats_dir=fullfile(subject.analysis_dir,subject.id,'stats');
end
if (~exist(subject.stats_dir,'dir'))
 mkdir(subject.stats_dir)
end

if (~isfield(subject,'sensor_stats'))
 subject.sensor_stats=fullfile(subject.stats_dir,['sensor_stats_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
end


%% Deal with trials selection
if (isfield(params,'nTrials'))
 nTrials=params.nTrials;
else
 nTrials=[]; % take all that are good
end

% check if we have a trial selection already
if (~isfield(subject,'SelectedTrials_file'))
 subject.SelectedTrials_file=fullfile(subject.analysis_dir,subject.id,'processed',['selected_trials_' subject.id '_' subject.exam_id '.mat']);
end
if (exist(subject.SelectedTrials_file,'file'))
 % if we have a file, load this, and be done
 disp('TrialSelection file found - loading')
 tsubj=load(subject.SelectedTrials_file,'subject');
 subject.SelectedTrials=tsubj.subject.SelectedTrials;
else
 % determine trials now
 badTrials = [];
 goodTrials = [1:length(data.trial)]; % start with all good
 
 
 if (isfield(params,'rejectEvents') && params.rejectEvents == 1)
  % check in dws_filt dataset not to miss atypically run processing
  if ~isfield(subject,'evt_timings_seconds') && ~isfield(subject,'evt_markerFile_notFound')
   subject_clean=load(subject.dws_filt_dataset,'subject');
   items=fieldnames(subject_clean.subject);
   for i=1:length(items)
    if length(items{i})>4 && strcmp(items{i}(1:4),'evt_')
     subject.(items{i})=subject_clean.subject.(items{i});
    end
   end
   if isfield(subject_clean.subject,'stamps') && isfield(subject_clean.subject.stamps,'readingevents')
    subject.stamps.readingevents=subject_clean.subject.stamps.readingevents;
   end
   clear subject_clean
  end
 end
 
  % call the central selection function now
 [ goodTrials, badTrials ] = nmri_trial_selector(subject,data,params);
  
 fprintf('\nTotal: GoodTrials, N=%d / BadTrials, N=%d\n',length(goodTrials),length(badTrials))
 subject.evt_goodTrials=goodTrials;
 subject.evt_badTrials=badTrials;
 
 if ~isempty(nTrials)
  % check if sufficient
  if (length(subject.evt_goodTrials)<nTrials)
   error(sprintf('Fewer good trials (%d) [after excluding events] in dataset than given in nTrials(%d)',length(subject.evt_goodTrials),nTrials))
  end
  % draw only from good trials then - make sure we are really random, use a
  % new rand stream
  s = RandStream('mt19937ar','Seed', seconds(round(milliseconds(second(datetime))*1000000)));
  subject.SelectedTrials = sort(datasample(s,subject.evt_goodTrials,nTrials,'Replace',false)); 
  fprintf('\nHave drawn N=%d random trials from all good trials now\n',length(subject.SelectedTrials))
 else
  % just take all good
  subject.SelectedTrials = goodTrials;
  fprintf('\nSelecting all good trials (N=%d) now\n',length(subject.SelectedTrials))
 end
  
 % and safe selection
 save(subject.SelectedTrials_file,'subject')
end

% we do not need trial markings any more now, remove to avoid Fieldtrip
% warnings
if isfield(data,'trial_markings')
 data=rmfield(data,'trial_markings');
end
if isfield(data,'trial_markings_sampleinfo')
 data=rmfield(data,'trial_markings_sampleinfo');
end


%% select just the wanted trials
cfg           = [];
cfg.trials    = subject.SelectedTrials;

% and also the non-bad channels - as of 01102019
if isfield(data,'bad_channels')
 good_channels={};
 for i=1:length(data.label)
  if ~any(strcmp(data.label{i},data.bad_channels))
   good_channels(end+1)=data.label(i);
  end
 end
 cfg.channel=good_channels;
 % we do not need bad_channels any more now, remove to avoid Fieldtrip
 % warnings
 if isfield(data,'bad_channels')
  data=rmfield(data,'bad_channels');
 end
 
end
data=ft_selectdata(cfg,data);


% Now demean finally and re-reference if EEG
if (strcmp(subject.dtype,'EEG'))
 cfg          = [];
 cfg.demean   = params.preproc_cfg.demean;
 cfg.reref      = 'yes';
 cfg.refchannel = 'all';
 data        = ft_preprocessing(cfg,data);
end



%% deal with layout mappping, we may have lost some channels along the way...

% get our layout
cfg = [];
cfg.layout=subject.layout;
cfg.skipscale='yes';
cfg.skipcomnt='yes';
layout=ft_prepare_layout(cfg);
layout.label=upper(layout.label); % make upper case

% get intersect
use_channels=intersect(layout.label,data.label);

% make a mapping to full layout

% now map from real labels to layout
layout_mapping=nan(length(use_channels),1);
for lx=1:length(use_channels)
 idx=find(strcmpi(use_channels{lx},layout.label));
 if length(idx)==1
  layout_mapping(lx)=idx;
 else
  idx
  error(['Could not find a unique mapping from channel (' use_channels{lx} ') in the data and the layout. This should not happen.'])
 end
end

% and remove all other channels
cfg           = [];
cfg.channel=use_channels;
data=ft_selectdata(cfg,data);


%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end


%% Now do full Fourier analysis 
% preallocate the freq structure
iFreq = 1;
fprintf('\nFourier transform - frequency: %s\n',params.freqsNames{iFreq}) 
cfg           = [];
cfg.method    = 'mtmfft';
cfg.output    = 'fourier';
cfg.pad       = 'nextpow2'; % recommended for speed
cfg.foi       = params.freqs(iFreq);
%cfg.trials    = subject.SelectedTrials;
cfg.tapsmofrq = params.tapsmofrq(iFreq); 
tmp   = ft_freqanalysis(cfg, data);
sensor_fourier = repmat(tmp,1,length(params.freqs));
    
for iFreq = 2:length(params.freqs)
 fprintf('\nFourier transform - frequency: %s\n',params.freqsNames{iFreq}) 
 cfg.foi       = params.freqs(iFreq);
 cfg.tapsmofrq = params.tapsmofrq(iFreq); 
 sensor_fourier(iFreq)   = ft_freqanalysis(cfg, data);
end
   
   
%% perform source space power


% calculate power and CSD in sensor space
for ff = 1:length(params.freqs)     
 fprintf('\ncalculate power and CSD - frequency: %s\n',params.freqsNames{iFreq}) 
 cfg           = [];
 cfg.method    = 'mtmfft';
 cfg.output    = 'powandcsd';
 cfg.pad       = 'nextpow2'; % recommended for speed
 cfg.foi       = sensor_fourier(ff).freq;
% cfg.trials    = subject.SelectedTrials;
 cfg.tapsmofrq = params.tapsmofrq(ff); 
 pow_and_csd_freq(ff)       = ft_freqanalysis(cfg, data);
end

% save
for ff = 1:length(params.freqs)
 % now remap according to layout
 metric=nan(length(layout.label),1);
 metric(layout_mapping,1)=pow_and_csd_freq(ff).powspctrm(:,1);
 save(fullfile(subject.stats_dir,['sensor_power_' subject.id '_' subject.exam_id '_' params.freqsNames{ff} '.mat' ]),'metric','layout')
end

save(subject.sensor_stats,'sensor_fourier','pow_and_csd_freq');

% make plot
disp('Now doing sensor space power plots...')   
allf=length(params.freqs);
cols=3;
rows=ceil(allf/cols);
hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','off'); 
if (~isfield(subject,'layout'))
 if (isfield(params,'layout'))
  subject.layout=params.layout;
 else
  error('no layout scheme specified in either subject or param')
 end
end
cfg=[];
cfg.gridscale = 200; %more beautiful
cfg.layout    = layout;
cfg.comment   = 'xlim';
for i=1:allf
subplot(rows,cols,i)
 title(params.freqsNames{i},'FontSize',12,'FontWeight','bold')
 ft_topoplotTFR(cfg,pow_and_csd_freq(i))
end
saveas(hFig,fullfile(subject.QCdir,['sensor_power_' subject.id '_' subject.exam_id '.png']),'png'); 
delete(hFig)



%% calculating connectivity metric by requested params


for mm=1:length(params.con_method)
 method = params.con_method{mm}; 
 cfg = [];
 cfg.method = method;
 if (strcmp(method, 'coh')) 
  cfg.complex = 'complex';
 end
 fprintf('\nDoing connectivity analysis for: %s\n', method)
 all_metric={};
 for ff = 1:length(params.freqs)
  fprintf('\nFrequency: %s\n', params.freqsNames{ff})
  tmp = ft_connectivityanalysis(cfg,sensor_fourier(ff));
  if (strcmp(method, 'wpli_debiased'))
   metric = tmp.wpli_debiasedspctrm;  
  elseif (strcmp(method, 'coh'))
   metric = tmp.cohspctrm;   
  end
  % map to layout
  all_metric{ff}=nan(length(layout.label));
  all_metric{ff}(layout_mapping,layout_mapping)=metric;
 end
 
 % save metric
 if (isfield(cfg,'complex') && strcmp(cfg.complex,'complex') && ~isreal(all_metric{ff}))
  % safe real  
  fprintf('Saving metric (real, %s)...\n', method)
  for ff = 1:length(params.freqs)
   metric=real(all_metric{ff});
   metric(eye(size(metric))==1)=NaN; % NaN on the identity line
   save(fullfile(subject.stats_dir,['sensor_' method '_real_' subject.id '_' subject.exam_id '_' params.freqsNames{ff} '.mat' ]),'metric','layout')
  end
  % Plot metric @sensors
  hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','Off'); 
  cfg=[];
  cfg.gridscale = 200; %more beautiful
  cfg.layout    = layout;
  cfg.comment   = 'xlim';
  
  cfg.colorbar  = 'EastOutside' ;
  % spoof in
  data_struct=pow_and_csd_freq;
  for i=1:allf
   metric=abs(real(all_metric{i}));
   data_struct(i).powspctrm=nanmean(metric(layout_mapping,layout_mapping))';
   subplot(rows,cols,i)
   title([strrep(method,'_','-') ' (real/abs) - ' params.freqsNames{i}],'FontSize',12,'FontWeight','bold')
   ft_topoplotTFR(cfg,data_struct(i))
  end
  saveas(hFig,fullfile(subject.QCdir,['sensor_con_' method '_real_' subject.id '_' subject.exam_id '.png']),'png'); 
  delete(hFig)


  % safe imaginary  
  fprintf('Saving metric (imaginary, %s)...\n', method)
  for ff = 1:length(params.freqs)
   metric=imag(all_metric{ff});
   %CAVE: imag of NaN is 0..., keep the NaN/non-available channels
   metric(isnan(all_metric{ff}))=NaN;
   metric(eye(size(metric))==1)=NaN; % NaN on the identity line
   save(fullfile(subject.stats_dir,['sensor_' method '_img_' subject.id '_' subject.exam_id '_' params.freqsNames{ff} '.mat' ]),'metric','layout')
  end
  % Plot metric @sensors
  hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','Off'); 
  cfg=[];
  cfg.gridscale = 200; %more beautiful
  cfg.layout    = layout;
  cfg.comment   = 'xlim';
  cfg.colorbar  = 'EastOutside' ;
  % spoof in
  data_struct=pow_and_csd_freq;
  for i=1:allf
   metric=abs(imag(all_metric{i}));
   data_struct(i).powspctrm=nanmean(metric(layout_mapping,layout_mapping))';
   subplot(rows,cols,i)
   title([strrep(method,'_','-') ' (img/abs) - ' params.freqsNames{i}],'FontSize',12,'FontWeight','bold')
   ft_topoplotTFR(cfg,data_struct(i))
  end
  saveas(hFig,fullfile(subject.QCdir,['sensor_con_' method '_img_' subject.id '_' subject.exam_id '.png']),'png'); 
  delete(hFig)

 else
  % safe as is
  fprintf('Saving metric (%s)...\n', method)
  for ff = 1:length(params.freqs)
   metric=all_metric{ff};
   metric(eye(size(metric))==1)=NaN; % NaN on the identity line
   save(fullfile(subject.stats_dir,['sensor_' method '_' subject.id '_' subject.exam_id '_' params.freqsNames{ff} '.mat' ]),'metric','layout')
  end
  % Plot metric @sensors
  hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','Off'); 
  cfg=[];
  cfg.gridscale = 200; %more beautiful
  cfg.layout    = layout;
  cfg.comment   = 'xlim';
  cfg.colorbar  = 'EastOutside' ;
  % spoof in
  data_struct=pow_and_csd_freq;
  for i=1:allf
   metric=all_metric{i};
   data_struct(i).powspctrm=nanmean(metric(layout_mapping,layout_mapping))';
   subplot(rows,cols,i)
   title([strrep(method,'_','-') ' - ' params.freqsNames{i}],'FontSize',12,'FontWeight','bold')
   ft_topoplotTFR(cfg,data_struct(i))
  end
  saveas(hFig,fullfile(subject.QCdir,['sensor_con_' method '_' subject.id '_' subject.exam_id '.png']),'png'); 
  delete(hFig)
 end
end

%% if we got this far, place a stamp for completed processing
subject=nmri_stamp_subject(subject,'processing_sensor',params);

end % final end

