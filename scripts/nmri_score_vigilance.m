function [ subject ] = nmri_score_vigilance(subject, data, params)
%[ subject ] = nmri_score_vigilance(subject, data, params)
%  
% This function as a wrapper to visually score vigilance and technical
% quality of EEG/MEG via the eeg_score viewer
% it requieres user interaction / presence
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file

% written by NF  03/2017


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

if (~exist('data','var') ) 
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

%% deal with selecting only the good trials, it can be a pain to rescore/view many bad trials and there is no real point anyways

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

data_r=data; % data_r is the same as data, but without the trial markings
if isfield(data_r,'trial_markings')
 data_r=removefields(data_r,'trial_markings');
end
if isfield(data_r,'trial_markings_sampleinfo')
 data_r=removefields(data_r,'trial_markings_sampleinfo');
end
if isfield(data_r,'bad_channels')
 data_r=removefields(data_r,'bad_channels');
end


cfg           = [];
cfg.trials    = goodTrials;
data_r=ft_selectdata(cfg,data_r);


% now re-generate trial markings
if isfield(data,'trial_markings')
 data_r.trial_markings=data.trial_markings(goodTrials,:);
 data_r.trial_markings_sampleinfo=data.trial_markings_sampleinfo(goodTrials,:);
end
if isfield(data,'bad_channels')
 data_r.bad_channels=data.bad_channels;
end


%% Now do the actual scoreing
h=msgbox({'Score each trial for vigilance.','Use keyboard (0-3) or from the dropdown.','You can also flag trials as technically bad,','these trials will be discarded from the final analysis as well',...
 '','Events are  also shown and can be checked / modified','',['Note: Minimum of trials needed = ' num2str(params.nTrials) ', currently = ' num2str(length(data_r.trial))],'Close the window when done'},'Trial inspection');
uiwait(h)

cfg = [];
cfg.score_vigilance=true;
cfg.score_technical=true;
cfg.allow_events=true; % do not allow to place events
cfg.select_montage=1; % take first of auto-detected
[data_r, subject]=eeg_score(cfg,data_r,subject);


% now copy back
data.trial_markings(goodTrials,:)=data_r.trial_markings;
if isfield(data_r,'bad_channels') 
 data.bad_channels=data_r.bad_channels;
end

%% check if todo=0, then stamp
if (sum(cellfun(@isempty,data_r.trial_markings(:,1)))==0)
 subject=nmri_stamp_subject(subject,'vigilance',params);
end

clear data_r

%% now place the scoring in the cleaned file(s)
if useICA
 nmri_write_dataset(subject.cleanICA_dataset,data,subject);
 % also place in clean
 data_clean=load(subject.clean_dataset,'data');
 % check length (just in case)
 if length(data.trial) ~= length(data_clean.data.trial)
  error('Not the same trial count between clean and cleanICA dataset, should not happen')
 end
 data_clean.data.trial_markings=data.trial_markings;
 data_clean.data.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
 if isfield(data,'bad_channels') 
  data_clean.data.bad_channels=data.bad_channels;
 end
 nmri_write_dataset(subject.clean_dataset,data_clean.data,subject);
 clear data_clean
else
 nmri_write_dataset(subject.clean_dataset,data,subject);
end

%% Check if events got changed - then put in dws_filt for reference
subject_dws=load(subject.dws_filt_dataset,'subject');
if ~isequal(subject_dws.subject,subject)
 fprintf('Updating subject events in %s\n',subject.dws_filt_dataset)
 save(subject.dws_filt_dataset,'subject','-append');
end

%% Check in hdm_lead
if isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')
 subject_hdm=load(subject.hdm_lead,'subject');
 if isfield(data,'trial_markings')
  subject_hdm.subject.data_info.trial_markings=data.trial_markings;
 end
 if isfield(data,'trial_markings_sampleinfo')
  subject_hdm.subject.data_info.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
 end
 if ~isequal(subject_hdm.subject,subject)
  update_markings(subject_hdm.subject,subject,subject.hdm_lead)
 end
end

%% Check in selected_trials
if isfield(subject,'SelectedTrials_file') && exist(subject.SelectedTrials_file,'file')
 subject_hdm=load(subject.SelectedTrials_file,'subject');
 if isfield(data,'trial_markings')
  subject_hdm.subject.data_info.trial_markings=data.trial_markings;
 end
 if isfield(data,'trial_markings_sampleinfo')
  subject_hdm.subject.data_info.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
 end
 if ~isequal(subject_hdm.subject,subject)
  update_markings(subject_hdm.subject,subject,subject.SelectedTrials_file)
 end
end



end


function update_markings(subject,subject_updated,fname)
 % update the relevant fields, leave rest intact
 items=fieldnames(subject_updated);
 for i=1:length(items)
  if length(items{i})>=4 && (strcmp(items{i}(1:4),'evt_'))
   subject.(items{i})=subject_updated.(items{i});
  end
  if length(items{i})>=15 && (strcmp(items{i}(1:15),'trial_markings'))
   subject.(items{i})=subject_updated.(items{i});
  end
 end
 fprintf('Updating events and trial markings in %s\n',fname)
 save(fname,'subject','-append');
end
