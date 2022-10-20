function[ subject, data] = cs_automatic_rejection(subject, data)
% This function runs automatic artefactrejection for magnetometers and gradiometers sepearately.
% Steps are based on https://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection/
% including identification of trials with muscle artefacts using z-values, additional lowpass filtering,
% and rejection of contaminated trials (detected in both channel types).
% Two datasets are saved in the final step: one with all original trials and indications which trials are bad (clean), 
% and one for which bad trials have been rejected directly (clean_rej). 
% This is needed for further steps like ICA-selection etc. Both datasets will contain both channel type data.
% Christina Stier, written and modified in 2020/22

%% load data and subject if needed
% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

% load data if needed
if (~isfield(subject,'dws_filt_dataset') || ~exist(subject.dws_filt_dataset,'file'))
 subject.dws_filt_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_' subject.id '_' subject.exam_id '.mat']);
 load(subject.dws_filt_dataset,'data');
else
 disp('Loading dws_filtered dataset...')
 load(subject.dws_filt_dataset,'data');
end


%% prepare the cleaning
bad_trials=zeros(1,length(data.trial));
bad_channels=zeros(length(data.label),1);

%% muscles MAG
% muscles magnetometers
cfg            = [];
cfg.trl        = data.trial;
%cfg.datafile   = 'ArtifactMEG.ds';
%cfg.headerfile = 'ArtifactMEG.ds';
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel      = 'MEG*1';
cfg.artfctdef.zvalue.cutoff       = 14;
cfg.artfctdef.zvalue.trlpadding   = 0;
cfg.artfctdef.zvalue.fltpadding   = 0;
cfg.artfctdef.zvalue.artpadding   = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter     = 'yes';
cfg.artfctdef.zvalue.bpfreq       = [110 140];
cfg.artfctdef.zvalue.bpfiltord    = 9;
cfg.artfctdef.zvalue.bpfilttype   = 'but';
cfg.artfctdef.zvalue.hilbert      = 'yes';
cfg.artfctdef.zvalue.boxcar       = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'no';

cfg.artfctdef.reject = 'complete'; 

[cfg, artifact_muscle_mag] = ft_artifact_zvalue(cfg,data);

%% save artifact_muscle mag
if (~isfield(subject,'artifact_muscle_mag'))
 subject.artifact_muscle_mag=fullfile(subject.analysis_dir,subject.id,'processed',['artifact_muscle_mag_' subject.id '_' subject.exam_id '.mat']);
end

disp('Writing artifact muscle mag file...')
save(subject.artifact_muscle_mag,'artifact_muscle_mag','subject')

%% muscles GRAD
% muscles gradiometers
cfg            = [];
cfg.trl        = data.trial;
%cfg.datafile   = 'ArtifactMEG.ds';
%cfg.headerfile = 'ArtifactMEG.ds';
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel      = {'MEG*2', 'MEG*3'};
cfg.artfctdef.zvalue.cutoff       = 14;
cfg.artfctdef.zvalue.trlpadding   = 0;
cfg.artfctdef.zvalue.fltpadding   = 0;
cfg.artfctdef.zvalue.artpadding   = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter     = 'yes';
cfg.artfctdef.zvalue.bpfreq       = [110 140];
cfg.artfctdef.zvalue.bpfiltord    = 9;
cfg.artfctdef.zvalue.bpfilttype   = 'but';
cfg.artfctdef.zvalue.hilbert      = 'yes';
cfg.artfctdef.zvalue.boxcar       = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'no';

cfg.artfctdef.reject = 'complete'; 

[cfg, artifact_muscle_grad] = ft_artifact_zvalue(cfg,data);

%% save artifact_muscle grad
if (~isfield(subject,'artifact_muscle_grad'))
 subject.artifact_muscle_grad=fullfile(subject.analysis_dir,subject.id,'processed',['artifact_muscle_grad_' subject.id '_' subject.exam_id '.mat']);
end

disp('Writing artifact muscle grad file...')
save(subject.artifact_muscle_grad,'artifact_muscle_grad','subject')
        

%% combine artifacts from mag & grad and save respective trials
artifact_muscle = [artifact_muscle_mag;artifact_muscle_grad];

if (~isfield(subject,'artifact_muscle_all'))
 subject.artifact_muscle_all=fullfile(subject.analysis_dir,subject.id,'processed',['artifact_muscle_all_' subject.id '_' subject.exam_id '.mat']);
end

disp('Writing artifact muscle all file...')
save(subject.artifact_muscle_all,'artifact_muscle','subject');

%% Additionally lowpass filter data (line noise at 50/150/200Hz!)
cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 70;
cfg.lpfiltord  = 1; % first order Butterworth
data = ft_preprocessing(cfg, data);

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
 % make sure we have enough columns (backward compatabilty)
 for i=2:4
  if size(data.trial_markings,2)<i
   data.trial_markings=[data.trial_markings cell(size(data.trial_markings,1),1)];
  end
 end
end

%% reject from data.trial directly
rcfg=[];
rcfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
%cfg.artfctdef.eog.artifact = artifact_EOG; %
%cfg.artfctdef.jump.artifact = artifact_jump;
rcfg.artfctdef.muscle.artifact = artifact_muscle;
data_no_artifacts = ft_rejectartifact(rcfg,data);


%% get trial index of bad trials and save in data.trialmarkings
idx = ~ismember(data.sampleinfo(:,1), data_no_artifacts.sampleinfo(:,1));
bad_trials = double(idx');

for i=1:length(data.trial)
 if bad_trials(i)
  % these trials have been rejected - mark as bad (false)
  data.trial_markings{i,2}=false;  
 end
end

%% make our output path and dir
if (~isfield(subject,'clean_dataset'))
 subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
end
if (~exist(fullfile(subject.analysis_dir,subject.id,'processed'),'dir'))
 mkdir(fullfile(subject.analysis_dir,subject.id,'processed'))
end

 %% save data with trial markings- prior to ICA
nmri_write_dataset(subject.clean_dataset,data,subject);
 
%% save clean data (trials rejected)
if (~isfield(subject,'clean_rej_dataset'))
subject.clean_rej_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_rej_' subject.id '_' subject.exam_id '.mat']);
end

data = data_no_artifacts;
nmri_write_dataset(subject.clean_rej_dataset,data,subject);
 
%% if we got this far, place a stamp for completed artifact rejection
subject=nmri_stamp_subject(subject,'artifactrejection',params);



