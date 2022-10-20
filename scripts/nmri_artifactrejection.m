function [ subject, data ] = nmri_artifactrejection(subject, data, params)
%[ subject, data ] =  nmri_artifactrejection(subject, data, params)
%  
% This function will do visual artifact on one subject
% it requieres user interaction / presenece
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 11/2016 
% updated to own viewers 03/2017


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

% make our output path and dir
if (~isfield(subject,'clean_dataset'))
 subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
end
if (~exist(fullfile(subject.analysis_dir,subject.id,'processed'),'dir'))
 mkdir(fullfile(subject.analysis_dir,subject.id,'processed'))
end

% check the cleaned version
if (~exist('data','var') ) 

if (isfield(params,'useICA_clean') && params.useICA_clean==1 && isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
 % have ICA
 if  isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
  %both 
  qcfg=[];
  qcfg.question={'Both a previously ICA-cleaned dataset and cleaned dataset were found.','',...
   'Do you want to work further on the ICA dataset (CleanedICA)',...
   'work with the cleaned dataset BEFORE ICA (Clean)',...
   'or start again with the raw/downsampled-fitlered dataset (Dws-filt)?'};
  qcfg.title='Reuse dataset?';
  qcfg.options={'CleanedICA','Cleaned','Dws-filt'};
  qcfg.mode='buttons';
  ref_dataset=nf_uidialog(qcfg);
 else
  % only ICA
  qcfg=[];
  qcfg.question={'A previously ICA-cleaned dataset was found.','',...
   'Do you want to work further on the ICA dataset (CleanedICA)',...
   'or start again with the raw/downsampled-fitlered dataset (Dws-filt)?'};
  qcfg.title='Reuse dataset?';
  qcfg.options={'CleanedICA','Dws-filt'};
  qcfg.mode='buttons';
  ref_dataset=nf_uidialog(qcfg);
 end
else
 % check for cleaned only, no ICA
 if  isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
  qcfg=[];
  qcfg.question={'A previously cleaned dataset was found.','',...
   'Do you want to work further on the this dataset (Cleaned)',...
   'or start again with the raw/downsampled-fitlered dataset (Dws-filt)?'};
  qcfg.title='Reuse dataset?';
  qcfg.options={'Cleaned','Dws-filt'};
  qcfg.mode='buttons';
  ref_dataset=nf_uidialog(qcfg);
 else
  ref_dataset='Dws-filt';
 end
end

if strcmpi(ref_dataset,'CleanedICA')
 disp('Loading previously ICA-cleaned dataset...')
 load(subject.cleanICA_dataset,'data');
 ref_data=subject.cleanICA_dataset;
elseif strcmpi(ref_dataset,'Cleaned')
 disp('Loading previously cleaned dataset...')
 load(subject.clean_dataset,'data');
 ref_data=subject.clean_dataset;
elseif strcmpi(ref_dataset,'Abort')
 disp('Stopping on user request...')
 return
else 
  % check if we have dws_filt_dataset (this is the minimum)
 if (~isfield(subject,'dws_filt_dataset') || ~exist(subject.dws_filt_dataset,'file'))
  error('Filtered and downsampled dataset not specified, run nmri_preproc first')
 else
  disp('Loading raw dataset...')
  load(subject.dws_filt_dataset,'data');
  ref_data=subject.dws_filt_dataset;
 end
end

end





% want ICA cleaning?
if (isfield(params,'useICA_clean') && params.useICA_clean==1)
 if (~isfield(subject,'cleanICA_dataset'))
  subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
 end
 if (~isfield(subject,'ICA_components'))
  subject.ICA_components=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_comp_' subject.id '_' subject.exam_id '.mat']);
 end
end

%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end



 %% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% make sure we have trial markings
if ~isfield(data,'trial_markings')
 data.trial_markings=cell(length(data.trial),4);
 % col1: sleep
 % col2: technical
 % col3: event
 % col4: rest/stimulation
 
 data.trial_markings_sampleinfo=cell(length(data.trial),2);
 % col1: sampleinfo start/stop
 % col2: seconds start/stop
 for i=1:length(data.trial)
  data.trial_markings_sampleinfo{i,1}=data.sampleinfo(i,:);
  data.trial_markings_sampleinfo{i,2}=data.sampleinfo(i,:)/data.fsample;
 end
end

if size(data.trial_markings,2)<4
 % make sure we have enough columns (backward compatability)
 for i=2:4
  if size(data.trial_markings,2)<i
   data.trial_markings=[data.trial_markings cell(size(data.trial_markings,1),1)];
  end
 end
end

 
 %% Find bad channels and bad trials using visual inspection - 1st pass
 % use the var option to discard noisy trials in the first round of
 % selection, also very important to check Kurtosis to reject bad
 % channels, unusually small Kurtosis indicates dead channels
 
 
 

% temporarily reject technically bad trials
rej_params=params;
if isfield(rej_params,'rejectEvents')
 % do not reject events now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectEvents');
end
if isfield(rej_params,'rejectVigilance')
 % do not reject vigilance now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectVigilance');
end

[ goodTrials, ~ ] = nmri_trial_selector(subject,data,rej_params);

data_r=data;
tcfg           = [];
tcfg.trials    = goodTrials;
if isfield(data,'bad_channels')
 good_channels={};
 for i=1:length(data.label)
  if ~any(strcmp(data.label{i},data.bad_channels))
   good_channels(end+1)=data.label(i);
  end
 end
 tcfg.channel=good_channels;
end



% remove trial markings to avoid Fieldtrip warnings

if isfield(data_r,'trial_markings')
 data_r=removefields(data_r,'trial_markings');
end
if isfield(data_r,'trial_markings_sampleinfo')
 data_r=removefields(data_r,'trial_markings_sampleinfo');
end
if isfield(data_r,'bad_channels')
 data_r=removefields(data_r,'bad_channels');
end
data_r=ft_selectdata(tcfg,data_r);

h=msgbox({'Check var and kurtosis (high in noisy trials/channels, low in dead channels), 1st pass','Start by rejecting trials, then channels (only if unavoidable)',['Note: Minimum of trials needed = ' num2str(params.nTrials) ', currently = ' num2str(length(data_r.trial))],'Press quit when done!! Do NOT hit the close window ''X'''},'Visual rejection');
uiwait(h)

ocfg          = [];
ocfg.method   = 'summary'; 
if (isfield(subject,'layout'))
 ocfg.layout  = subject.layout ;
elseif (strcmp(subject.dtype,'MEG'))
 % preset for MEG
 ocfg.layout    = 'CTF275.lay';
else
 error('Layout could not be found...')
end

data_r         = ft_rejectvisual(ocfg,data_r); % press quit when done



% as of 05/2019 changed to no really delete the data/trials, but rather to
% mark as technically bad

for i=1:length(data.trial)
 found_t=find(data.sampleinfo(i,1)==data_r.sampleinfo(:,1));
 if isempty(found_t)
  % this trials seems to have been rejected - mark as bad (false)
  data.trial_markings{i,2}=false;  
 end
end


% check if channels were rejected
bad_channels={};
for i=1:length(data.label)
 if ~any(strcmp(data.label{i},data_r.label))
  bad_channels(end+1)=data.label(i);
 end
end

% add to bad channels in data
data.bad_channels=bad_channels;



% we do not do a real rejection any more since 01/10/2019
% % rejected channels will be removed from the data
% tcfg          = [];
% tcfg.channel   = data_r.label; 
% 
% % remove trial markings to avoid Fieldtrip warnings
% data_r=data;
% if isfield(data_r,'trial_markings')
%  data_r=removefields(data_r,'trial_markings');
% end
% if isfield(data_r,'trial_markings_sampleinfo')
%  data_r=removefields(data_r,'trial_markings_sampleinfo');
% end
% 
% data_r=ft_selectdata(tcfg,data_r);
% 
% if isfield(data,'trial_markings')
%  data_r.trial_markings=data.trial_markings;
%  data_r.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
% end
% 
% data=data_r;

% temporarily reject technically bad trials
rej_params=params;
if isfield(rej_params,'rejectEvents')
 % do not reject events now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectEvents');
end
if isfield(rej_params,'rejectVigilance')
 % do not reject vigilance now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectVigilance');
end

[ goodTrials, ~ ] = nmri_trial_selector(subject,data,rej_params);

if (length(goodTrials)<params.nTrials)
 warning(sprintf('Fewer trials (%d) in dataset than requested in nTrials(%d)',length(goodTrials),params.nTrials))
end

fprintf('Channels: Good, N=%d / Bad, N=%d\n',length(data.label)-length(unique(data.bad_channels)),length(unique(data.bad_channels)))

%% now do visual inspection trial-by-trial

% deal with the trial marking fields again
data_r=data; % data_r is the same as data, but without the trial markings and bad channels
if isfield(data_r,'trial_markings')
 data_r=removefields(data_r,'trial_markings');
end
if isfield(data_r,'trial_markings_sampleinfo')
 data_r=removefields(data_r,'trial_markings_sampleinfo');
end
if isfield(data_r,'bad_channels')
 data_r=removefields(data_r,'bad_channels');
end


% make common AVG montage (temporary), if EEG... just for eeg_score reading
if (strcmp(subject.dtype,'EEG')) || strcmp(subject.dtype,'EEG_invasive')
 fprintf('Making a temporary common AVG...\n')
 cfg          = [];
 cfg.demean   = params.preproc_cfg.demean;
 cfg.reref      = 'yes';
 cfg.refchannel = 'all';
 data_avg        = ft_preprocessing(cfg,data_r);
else
 data_avg=data_r;
end

clear data_r

% now copy the trial_markings back
if isfield(data,'trial_markings')
 data_avg.trial_markings=data.trial_markings;
 data_avg.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
end
if isfield(data,'bad_channels')
 data_avg.bad_channels=data.bad_channels;
end

% changed to use our own data browser
h=msgbox({'Check each trial for viusal artifacts / movement.','Mark bad trials (backspace) or good trials (enter), or from the dropdown.',...
 'The whole trial will be discarded then.','Also check if stimuli are detected correctly','','Vigilance scoring is blocked here (should be done after the cleaning)',...
 'Events are shown, but it is not possible to place new events now','',...
 ['Note: Minimum of trials needed = ' num2str(params.nTrials) ', currently = ' num2str(length(goodTrials))],'Close the window when done'},'Trial inspection');
uiwait(h)

% now call the viewer
cfg = [];
cfg.score_vigilance=false;
cfg.score_technical=true;
cfg.allow_events=false; % do not allow to place events
cfg.select_montage=1; % take first of auto-detected
cfg.hide_bad_channels=0; % do not hide bad channels by default
cfg.plot_layout=1; % plot the layout as default
data_avg=eeg_score(cfg,data_avg,subject);

% save trial_markings back to data
if isfield(data_avg,'trial_markings')
 data.trial_markings=data_avg.trial_markings;
 data.trial_markings_sampleinfo=data_avg.trial_markings_sampleinfo;
end
if isfield(data_avg,'bad_channels')
 data.bad_channels=data_avg.bad_channels;
end


% now reject the trials marked as bad  (we do not have vigilance scored at this stage usually)
% so will use technical and stimuli

% temporarily reject technically bad trials
rej_params=params;
if isfield(rej_params,'rejectEvents')
 % do not reject events now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectEvents');
end
if isfield(rej_params,'rejectVigilance')
 % do not reject vigilance now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectVigilance');
end

[ goodTrials, ~ ] = nmri_trial_selector(subject,data,rej_params);

clear data_r data_avg
if (length(goodTrials)<params.nTrials)
 warning(sprintf('Fewer trials (%d) in dataset than given in nTrials(%d)',length(goodTrials),params.nTrials))
end


fprintf('Channels: Good, N=%d / Bad, N=%d\n',length(data.label)-length(unique(data.bad_channels)),length(unique(data.bad_channels)))



%% Now do the 2nd pass of summary stats
% display the summary again and save the resulting data
% you may discard further trials if the trial numbers are significantly

h=msgbox({'2nd pass of summary stats',['Note: Minimum of trials needed = ' num2str(params.nTrials) ', currently = ' num2str(length(goodTrials))],'Close window or press quit when done!! Do NOT hit the close window ''X'''},'Visual rejection');
uiwait(h)

% remove trial markings to avoid Fieldtrip warnings
data_r=data;
if isfield(data_r,'trial_markings')
 data_r=removefields(data_r,'trial_markings');
end
if isfield(data_r,'trial_markings_sampleinfo')
 data_r=removefields(data_r,'trial_markings_sampleinfo');
end
if isfield(data_r,'bad_channels')
 data_r=removefields(data_r,'bad_channels');
end



% temporarily reject technically bad trials
rej_params=params;
if isfield(rej_params,'rejectEvents')
 % do not reject events now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectEvents');
end
if isfield(rej_params,'rejectVigilance')
 % do not reject vigilance now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectVigilance');
end

[ goodTrials, ~ ] = nmri_trial_selector(subject,data,rej_params);

tcfg           = [];
tcfg.trials    = goodTrials;
if isfield(data,'bad_channels')
 good_channels={};
 for i=1:length(data.label)
  if ~any(strcmp(data.label{i},data.bad_channels))
   good_channels(end+1)=data.label(i);
  end
 end
 tcfg.channel=good_channels;
 fprintf('Channels: Good, N=%d / Bad, N=%d\n',length(good_channels),length(data.bad_channels))
end
data_r=ft_selectdata(tcfg,data_r);



 % make common AVG montage (temporary), if EEG...
if (strcmp(subject.dtype,'EEG')) || strcmp(subject.dtype,'EEG_invasive')
 fprintf('Making a temporary common AVG...\n')
 cfg          = [];
 cfg.demean   = params.preproc_cfg.demean;
 cfg.reref      = 'yes';
 cfg.refchannel = 'all';
 data_avg        = ft_preprocessing(cfg,data_r);
else
 data_avg=data_r;
end



data_r         = ft_rejectvisual(ocfg,data_avg); % press quit when done



% as of 05/2019 changed to no really delete the data/trials, but rather to
% mark as technically bad

for i=1:length(data.trial)
 found_t=find(data.sampleinfo(i,1)==data_r.sampleinfo(:,1));
 if isempty(found_t)
  % this trials seems to have been rejected - mark as bad (false)
  data.trial_markings{i,2}=false;  
 end
end


% check if channels were rejected
bad_channels={};
for i=1:length(data.label)
 if ~any(strcmp(data.label{i},data_r.label))
  bad_channels(end+1)=data.label(i);
 end
end

% add to bad channels in data
data.bad_channels=bad_channels;



% now make the final common AVG montage if EEG...
if (strcmp(subject.dtype,'EEG')) || strcmp(subject.dtype,'EEG_invasive')
 fprintf('Making the final common AVG...\n')
 
 % remove trial markings to avoid Fieldtrip warnings
 data_r=data;
 if isfield(data_r,'trial_markings')
  data_r=removefields(data_r,'trial_markings');
 end
 if isfield(data_r,'trial_markings_sampleinfo')
  data_r=removefields(data_r,'trial_markings_sampleinfo');
 end
 if isfield(data_r,'bad_channels')
  data_r=removefields(data_r,'bad_channels');
 end
 
 cfg          = [];
 cfg.demean   = params.preproc_cfg.demean;
 cfg.reref      = 'yes';
 if isfield(data,'bad_channels')
  good_channels={};
  for i=1:length(data.label)
   if ~any(strcmp(data.label{i},data.bad_channels))
    good_channels(end+1)=data.label(i);
   end
  end
  cfg.refchannel=good_channels;
 else
  cfg.refchannel = 'all';
 end
 data_r        = ft_preprocessing(cfg,data_r);
 
 % copy back the trial markings
 if isfield(data,'trial_markings')
  data_r.trial_markings=data.trial_markings;
  data_r.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
 end
 if isfield(data,'bad_channels')
  data_r.bad_channels=data.bad_channels;
 end
 data=data_r;
 clear data_r
end



% now make a final count
rej_params=params;
if isfield(rej_params,'rejectEvents')
 % do not reject events now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectEvents');
end
if isfield(rej_params,'rejectVigilance')
 % do not reject vigilance now, only stimuli and technical
 rej_params=rmfield(rej_params,'rejectVigilance');
end

[ goodTrials, ~ ] = nmri_trial_selector(subject,data,rej_params);

if (length(goodTrials)<params.nTrials)
 warning(sprintf('Fewer trials (%d) in dataset than requested in nTrials(%d)',length(goodTrials),params.nTrials))
end


fprintf('Channels: Good, N=%d / Bad, N=%d\n',length(data.label)-length(unique(data.bad_channels)),length(unique(data.bad_channels)))


%% check for change and prompt for overwrite
if exist('ref_dataset','var') && ~strcmp(ref_dataset,'Dws-filt') && exist('ref_data','var') && exist(ref_data,'file')
 % load the original and check for change
 old_data=load(ref_data);
 if ~isequal(data.label,old_data.data.label) || ~isequal(data.trial_markings,old_data.data.trial_markings) || ~isequal(size(data.trial),size(old_data.data.trial)) 
  % there was a change, prompt for overwrite

  if strcmp(ref_dataset,'Cleaned') && (isfield(params,'useICA_clean') && params.useICA_clean==1 && isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
   % changed the cleaned dataset and there is an ICA cleaned
   qcfg=[];
   qcfg.question={['There was a change in your dataset (' ref_dataset ')'],'',...
   'If you change this dataset, the ICA needs to be re-run',...
   'Do you really want to overwrite the dataset AND delete the ICA-cleaned version?'};
   qcfg.title='Confirm Overwrite';
   qcfg.options={'Yes','No'};
   qcfg.default={'No'};
   qcfg.mode='buttons';
   button=nf_uidialog(qcfg);
   if strcmpi(button,'Yes')
    disp('Deleting ICA...')
    delete(subject.cleanICA_dataset);
    if (isfield(subject,'ICA_components') && exist(subject.ICA_components,'file'))
     delete(subject.ICA_components);
    end
   end
   % then continue with the rest
  else
   % ICA, just reset this
   qcfg=[];
   qcfg.question={['There was a change in your dataset (' ref_dataset ')'],'',...
   'Do you really want to overwrite the dataset?'};
   qcfg.title='Confirm Overwrite';
   qcfg.options={'Yes','No'};
   qcfg.default={'No'};
   qcfg.mode='buttons';
   button=nf_uidialog(qcfg);
  end
  if ~strcmpi(button,'Yes')
   disp('Changes are not saved on user request...Exiting')
   return
  end
 else
  disp('There were no changes made (it seems at least)...Exiting')
  return  
 end
end


%% if we got this far, place a stamp for completed artifact rejection
subject=nmri_stamp_subject(subject,'artifactrejection',params);


 %% save cleaned data - prior to ICA
 if exist('ref_data','var') && ~strcmp(ref_dataset,'Dws-filt') 
  nmri_write_dataset(ref_data,data,subject);
 else
  nmri_write_dataset(subject.clean_dataset,data,subject);
 end






%% now do the ICA cleaning if this is wanted
if (isfield(params,'useICA_clean') && params.useICA_clean==1) && (~exist('ref_dataset','var') || strcmp(ref_dataset,'Dws-filt') || strcmp(ref_dataset,'Cleaned'))  % ask if ICA should be run
 qcfg=[];
 qcfg.question={'Do you want to start the ICA processing and review now?','',...
 'Can take some time, usually be more efficient to run from GUI'};
 qcfg.title='Run ICA correction';
 qcfg.options={'Yes','No'};
 qcfg.default={'No'};
 qcfg.mode='buttons';
 button=nf_uidialog(qcfg);
 
 if strcmpi(button,'Yes') 
  subject=nmri_artifactrejection_estimateICA(subject);
  subject=nmri_artifactrejection_reviewICA(subject);
 end
end


%% check if we want a QC power plot
if (isfield(params,'QC_plots') && params.QC_plots==1 && (~exist('ref_dataset','var') || strcmp(ref_dataset,'CleanedICA')))
 % make plot for ICA
 qcfg=[];
 qcfg.question={'Do you want to generate a QC power plot now?',...
 'Can take a few minutes, depending on the size of the dataset/channels'};
 qcfg.title='Generate QC plot?';
 qcfg.options={'Yes','No'};
 qcfg.default={'No'};
 qcfg.mode='buttons';
 button=nf_uidialog(qcfg);
 
 if strcmpi(button,'Yes')
  
  % call the central QC power plot script
  nmri_qcplot_power_sensor(subject,data,params,['pow_sens_plot_postICA_' datestr(now,'yyyymmdd')]);
  
 end
end

