function[subject] = cs_visual_datacheck(subject)
% This function loads the ICA-cleaned dataset, which can be viewd using
% eeg_score. Manual rejection of trials possible using backshift, which
% will be saved in trial markings. Finally, all additionally rejected
% trials will be saved in the cleanICA and clean_rej datasets

% Written and modified by Christina Stier, 2020/22

%% load cleaned dataset (all trials)
if (~isfield(subject,'cleanICA_dataset'))
 subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
end

load(subject.cleanICA_dataset);

%% also load rej-cleaned dataset 
if (~isfield(subject,'clean_rej_dataset'))
 subject.clean_rej_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_rej_' subject.id '_' subject.exam_id '.mat']);
end

d_rej = load(subject.clean_rej_dataset);
data_rej = d_rej.data;


%% add trial markings
% if ~isfield(data,'trial_markings')
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
% end

if size(data.trial_markings,2)<4
 % make sure we have enough columns (backward compatabilty)
 for i=2:4
  if size(data.trial_markings,2)<i
   data.trial_markings=[data.trial_markings cell(size(data.trial_markings,1),1)];
  end
 end
end

%% view data
% call the subject and params include
nmri_include_read_ps

% score
cfg = [];
cfg.score_vigilance=true;
cfg.score_technical=true;
cfg.allow_events=false; % do not allow to place events
cfg.select_montage=1; % take first of auto-detected
cfg.hide_bad_channels=0; % do not hide bad channels by default
cfg.plot_layout=1; % plot the layout as default
data_avg=eeg_score(cfg,data,subject);

% save trial_markings back to data
if isfield(data_avg,'trial_markings')
 data.trial_markings=data_avg.trial_markings;
 data.trial_markings_sampleinfo=data_avg.trial_markings_sampleinfo;
end
if isfield(data_avg,'bad_channels')
 data.bad_channels=data_avg.bad_channels;
end

% save trial_markings back to data rej
if isfield(data_avg,'trial_markings')
 data_rej.trial_markings=data_avg.trial_markings;
 data_rej.trial_markings_sampleinfo=data_avg.trial_markings_sampleinfo;
end
if isfield(data_avg,'bad_channels')
 data_rej.bad_channels=data_avg.bad_channels;
end

[ goodTrials, ~ ] = nmri_trial_selector(subject,data);

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
% data=ft_selectdata(tcfg,data); % reject in ICAclean dataset
data_rej=ft_selectdata(tcfg,data_rej); % reject in clean rej dataset

if (length(goodTrials)<params.nTrials)
 warning(sprintf('Fewer trials (%d) in dataset than given in nTrials(%d)',length(goodTrials),params.nTrials))
 dropout_file = fullfile(subject.analysis_dir, 'subjects_dropout.mat');
 load(dropout_file);
 reason = 'too few trials';
 rej_subjects{end+1,1} = subject.id;
 rej_subjects{end,2} = reason;
 save('subjects_dropout.mat', 'rej_subjects')
end

fprintf('Channels: Good, N=%d / Bad, N=%d\n',length(data.label)-length(unique(data.bad_channels)),length(unique(data.bad_channels)))

%% save changed datasets
nmri_write_dataset(subject.cleanICA_dataset,data,subject);

data = data_rej; % overwrite  data with dataset in which bad trials were rejected
% remove trial markings, which have changed while selecting only good
% trials
if isfield(data,'trial_markings')
 data=rmfield(data,'trial_markings');
end
if isfield(data,'trial_markings_sampleinfo')
 data=rmfield(data,'trial_markings_sampleinfo');
end

nmri_write_dataset(subject.clean_rej_dataset,data,subject);

%% Decide on new round of ICA
% load cleaned dataset and again check
button=questdlg({'Do you want to restart the ICA cleaning?'},'Start ICA?','No');
if strcmpi(button,'Yes')
 subject = cs_nmri_artifactrejection_estimateICA_grad(subject);
 subject = cs_nmri_artifactrejection_estimateICA_mag(subject);
 subject = cs_find_ecg_2nd(subject);
 subject = cs_find_eog_2nd(subject);
 subject = cs_nmri_artifactrejection_rejectICA(subject);
 
 fprintf('Subject %s has been fully preprocessed [second run]',subject.id)
end


end
