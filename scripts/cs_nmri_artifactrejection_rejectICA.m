function [ subject, data ] = cs_nmri_artifactrejection_rejectICA(subject, data, params)
%[ subject, data ] =  nmri_artifactrejection(subject, data, params)
%  
% This is a function to reject the ICA components based on previously
% selected ICA components (based on similarity with EOG and ECG signals).
% This time the cleaned data of subjects are used, for which bad trials
% were rejected (clean_rej).
% Output will be a cleanedICA-dataset per subject and QC plot of the
% ICA-cleaned data (power at Grads and Mags)
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 09/2019 and modified by Christina Stier 2020/22



% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps


% check if we have cleaned dataset (bad trials rejected)
subject.clean_rej_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_rej_' subject.id '_' subject.exam_id '.mat']);
load(subject.clean_rej_dataset,'data');

% load the MAG components
if (~isfield(subject,'ICA_components_mag'))
  subject.ICA_components_mag=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_comp_mag_' subject.id '_' subject.exam_id '.mat']);
end
 
disp('Load ICA comp mag file...')
mag = load(subject.ICA_components_mag, 'comp');

if (~isfield(subject,'ICA_components_grad'))
  subject.ICA_components_grad=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_comp_grad_' subject.id '_' subject.exam_id '.mat']);
end
 
disp('Load ICA comp grad file...')
grad= load(subject.ICA_components_grad, 'comp');

%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );

%% now do the ICA cleaning if this is wanted
% if (isfield(params,'useICA_clean') && params.useICA_clean==1)

%% Load selected components for GRAD and MAG based on correlation / coherence analysis 
% Loop over channel type and and component selections (ecg/eog) and create a
% variable comprising all components selected for each channel type
% (comp_selected.(channel type)

modality = {'ecg','eog'};
chan = {'mag', 'grad'};
comp_chan = {};

for c = 1:length(chan)
 for m = 1:length(modality)
  filename = fullfile(subject.analysis_dir,subject.id,'processed', [modality{m}, '_', chan{c}, '_comp_', subject.id '_' subject.exam_id '.mat']);
  value = load(filename);
  comp_chan{m} = unique(cell2mat(struct2cell(value))); % get rid of repeated values
 end
 comp_chan_all = horzcat(comp_chan{:});
 comp_selected.(chan{c}) = comp_chan_all;
end

%% now get data for selected components and save in a file for further plotting later on 
chan = {'mag', 'grad'};
comp_data = {};
comp_data{1} =  mag.comp; % get dataset ICA components MAG
comp_data{2} = grad.comp; % get dataset ICA components GRAD

for c = 1:length(chan) % loop over channel type
  
 label_comp = {}; % store labels of selected components
 
 for sc = 1:length(comp_selected.(chan{c}))
  label = comp_selected.(chan{c})(sc);
  label_comp{sc} = mag.comp.label{label}; % select from one comp dataset
 end
 
 cfg = [];
 cfg.channel = label_comp;
 sel_comp_data = ft_selectdata(cfg, comp_data{c}); % get data for those components
 
 datname = fullfile([chan{c} '_selected_comp_data']);
 subject.(datname)=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_sel_comp_' chan{c} '_' subject.id '_' subject.exam_id '.mat']);
 nmri_write_dataset(subject.(datname),sel_comp_data,subject); % write dataset with selected components only
end

%% reject components and backproject
% take only the good trial selection, which should be identical for MAGs and GRADs.
% call the subject and params include
nmri_include_read_ps
% call the central selection function
badTrials = [];
goodTrials = [1:length(data.trial)]; % start with all good
[goodTrials, badTrials] = nmri_trial_selector(subject,data,params);

% get MAG data
cfg = [];
cfg.channel = {'MEG*1'};
cfg.trials = goodTrials;
data_mag = ft_selectdata(cfg, data);

cfg = [];
cfg.component = comp_selected.mag;
data_mag = ft_rejectcomponent(cfg,mag.comp,data_mag);

%get GRAD data
cfg = [];
cfg.channel = {'MEG*2', 'MEG*3'};
cfg.trials = goodTrials;
data_grad = ft_selectdata(cfg, data);

cfg = [];
cfg.component = comp_selected.grad;
data_grad = ft_rejectcomponent(cfg,grad.comp,data_grad);

% append data
fulldata = ft_appenddata([],data_mag,data_grad);
fulldata.grad = data_mag.grad;
%eeg_score = ([],fulldata,subject)

 %% if we got this far, place a stamp for completed artifact rejection
 subject=nmri_stamp_subject(subject,'cleanICA',params);
  
 %% save ICA cleaned data - after ICA
 if (~isfield(subject,'cleanICA_dataset'))
 subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
end
 nmri_write_dataset(subject.cleanICA_dataset,fulldata,subject);
 
 %% check if we want a QC power plot
if (isfield(params,'QC_plots') && params.QC_plots==1)
%  % make plot
%  
%  button=questdlg({'Do you want to generate a QC power plot now?','Can take a few minutes'},'Generate QC plot?','No');
%  if strcmpi(button,'Yes')
%   
  % call the central QC power plot script
  cs_nmri_qcplot_power_mag(subject,data,params,['mag_sensor_powerplot_postICA_' datestr(now,'yyyymmdd')]);
  cs_nmri_qcplot_power_grad(subject,data,params,['grad_sensor_powerplot_postICA_' datestr(now,'yyyymmdd')]);
%  end
end
close all

