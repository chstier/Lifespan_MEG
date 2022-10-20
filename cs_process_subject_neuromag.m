%%% This script contains all single subject preprocessing steps
%%% of CamCAN-data (resting-state, maxfiltered) 

% written by Christina Stier, 2020/2022

%% Load full list of CamCAN subjects
subjlist = importdata('participants.tsv');
subjects = subjlist.textdata(2:end,1);

%% Selection of subjects 
% Either use all subjects available, make preselections of subjects 
% or randomly draw subjects (for age groups) 

% randomly draw n = 2 subjects for testing
rng(0,'twister');
r = randi([0 650],2,1);
selected_subjects = subjects(r,:); 

% orelse select specific age-groups and save
% save('selected_subjects_5.mat', 'selected_subjects')

% orelse use all subjects in list
% selected_subjects = subjects;

clear subjects

%% Preprocessing

for i = 1:length(selected_subjects)
 subject = selected_subjects(i);
 
 subj_id = char(selected_subjects(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 
 % indicate where data are stored 
 path_directory = '/home/uni10/nmri/projects/cstier/aging/cc700/meg/pipeline/release004/BIDS_20190411/meg_rest_mf/';
 
 % read data and prepare subject-struct and preprocess 
 if (~exist(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.m'),'file'))
  subject = cs_nmri_read_subject(subj_id, exam_id, path_directory);
 end
 
 % run preprocessing
 % a downsampled, filtered, and trial segmented dataset will be created for the subject (details are found in the manuscript)
 % change parameters for preprocessing in analysis_params.m if needed
 fprintf('Preprocessing subject %s', subj_id)
 
 dws_data = fullfile(analysis_dir, subj_id, 'processed', ['dws_filt_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat']);
 if exist(dws_data,'file')
  load(dws_data);
  subject = nmri_load_subject_most_advanced(subject);
 else
  subject = nmri_preproc(subject);
 end
 
 % run automatic artefactrejection (muscles), compute and select ICA components, which show similarity with ECG/EOG, 
 % and reject them accordingly
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['clean_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_automatic_rejection(subject);  
 end
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['ICA_comp_grad_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_nmri_artifactrejection_estimateICA_grad(subject); % estimate 50 components
 end
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['ICA_comp_mag_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_nmri_artifactrejection_estimateICA_mag(subject); % estimate 50 component
 end
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['dws_filt_physio_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_preproc_physio(subject);
 end
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['ecg_mag_comp_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_find_ecg(subject);
 end
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['eog_mag_comp_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_find_eog(subject);
 end
 if (~exist(fullfile(analysis_dir, subj_id, 'processed', ['cleanICA_' subj_id '_ses-rest_meg_ses-rest_task-rest_proc-sss.mat'])))
  subject = cs_nmri_artifactrejection_rejectICA(subject);
 end
 if (~exist(fullfile(subject.QCdir,['grad_sensor_freqplot_' subj_id '_' subject.exam_id '.png']),'file'))
  disp('Now doing frequency plots for QC') 
  subject = cs_freqplot(subject);
 end
 clear subject
end

%% Quality check of preprocessing
% these steps load the ICA components + time series and sensor power plots for a check, and later the MEG data 
% for further manual rejection of trials if needed. Changes will be saved in the ICAcleaned and clean_rej datasets.

for i = 1:length(selected_subjects)
 
 % load subject struct
 subject = selected_subjects(i);
 subj_id = char(selected_subjects(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 load(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.mat'))
 subject = nmri_load_subject_most_advanced(subject);
 
 % check rejected ICA components
 subject = cs_check_selected_components(subject); 
 %subject = cs_nmri_artifactrejection_reviewICA_grad(subject); % code for manual rejection if necessary
 %subject = cs_nmri_artifactrejection_reviewICA_mag(subject); % code for manual rejection if necessary
 %subject = cs_nmri_artifactrejection_rejectICA2(subject); % code for manual rejection if necessary
 
 % check power plots?
 button=questdlg({'Do you want to review QC sensor power plots now?'},'Review QC plot?','No');
 if strcmpi(button,'Yes')
  qc_sensors = dir(fullfile(subject.QCdir, '*sensor_*'));
  for p = 1:length(qc_sensors)
   image = qc_sensors(p).name;
   im_show = imread(fullfile(subject.QCdir, image));
   f = figure;
   imshow(im_show);
   uiwait(f)
  end
 end

 % load cleaned dataset and again check
 button=questdlg({'Do you want to review ICA-cleaned datasets now?'},'Review QC plot?','No');
 if strcmpi(button,'Yes') 
  subject = cs_visual_datacheck(subject);
 end
 
 fprintf('Quality check of %s finished',subject.id)
 
 button=questdlg({'Review done'},'Keep subject?','No');
 
 if strcmpi(button,'No')
   dropout_file = fullfile(subject.analysis_dir, 'subjects_dropout.mat');
   
   if exist(dropout_file)
    load(dropout_file);
    reason = 'too noisy';
    rej_subjects{end+1,1} = subject.id;
    rej_subjects{end,2} = reason;
   else
    reason = 'too noisy';
    rej_subjects = {};
    rej_subjects{end+1,1} = subject.id;
    rej_subjects{end,2} = reason;
   end
    save('subjects_dropout.mat', 'rej_subjects')
   
 end
 
  close all
end

%% Alignment and headmodel
% For the final processing, we only include subjects with at least 30 trials (300s) of data 
% Hence, discard all subjects who have less than 30 trials from list but also move full 
% subject-directory to a drop-out folder 

% Note: for the headmodel making, the Freesurfer-reconstruction of individual brains and surface-based mapping (SUMA, AFNI) should be ready
% (Saad, Z. S., & Reynolds, R. C. (2012). Suma. Neuroimage, 62(2), 768-773. https://doi.org/10.1016/j.neuroimage.2011.09.016) 
% Add SUMA processing to the standard Freesurfer procedure using a similar command as below (ld = density factor)
% @SUMA_Make_Spec_FS -sid $fs -fspath ${SUBJECTS_DIR}/${fs} -GIFTI -ld 10 
% More information can be found here: https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/fs/fs_fsprep.html#tut-fs-fsprep

analysis_dir = pwd;
dropout_file = fullfile(analysis_dir, 'subjects_dropout.mat');

if exist(dropout_file)
 load(dropout_file);
else
 rej_subjects = {};
 save('subjects_dropout.mat', 'rej_subjects')
 load(dropout_file);
end

if ~isempty(rej_subjects)
 good_subjects = setdiff(selected_subjects, rej_subjects(:,1)); 
 dropout_folder = '/home/uni10/nmri/projects/cstier/aging_processing/dropout';

 for s = 1:length(rej_subjects(:,1))
  if exist(char(rej_subjects(s,1)), 'file')
   if exist(dropout_folder)
    movefile(char(rej_subjects(s,1)), dropout_folder)
   else 
    mkdir dropout
    movefile(char(rej_subjects(s,1)), dropout_folder)
   end
  end
 end

else
 good_subjects = selected_subjects;
end

% Align anatomy, make headmodels and run source processing for all subjects
% with at least 30 clean trials

for i = 1:length(good_subjects)
 
 % load subject struct
 subject = char(good_subjects(i));
 subj_id = char(good_subjects(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 load(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.mat'))
 subject = nmri_load_subject_most_advanced(subject);
 
 % specify freesurfer path (note: the freesurfer reconstruction should be done at this point)
 fsfolder = '/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0';
 anatomy = fullfile(fsfolder, [ subject.id '_anat_T2_T1']);
 
 % specify folder where anatomical CamCAN-data and landmark files are stored
 anat_orig_folder = '/home/uni10/nmri/projects/cstier/aging_anat';
 
 if exist (anatomy, 'dir')
 
  % align anatomy
  subject = cs_nmri_align_mri(subject, anat_orig_folder); 
  
  % make headmodel
  subject = cs_nmri_make_hdm_suma(subject);

 else
  dropout_file = fullfile(subject.analysis_dir, 'subjects_dropout.mat');
   
  if exist(dropout_file)
    load(dropout_file);
  else
    rej_subjects = {};
    save('subjects_dropout.mat', 'rej_subjects')
    load(dropout_file);
  end
   reason = 'no anatomy';
   rej_subjects{end+1,1} = subject.id;
   rej_subjects{end,2} = reason;
   save('subjects_dropout.mat', 'rej_subjects')
  
 end

 clear subject
 end

%% Quality check headmodel / projection

for i = 1:length(good_subjects)
 
 % load subject struct
 subject = good_subjects(i);
 subj_id = char(good_subjects(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 load(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.mat'))
 subject = nmri_load_subject_most_advanced(subject);
 
 % get hdm plots and source plots from QC folder
 qc_hdm = dir(fullfile(subject.QCdir, 'hdm_*.fig'));
 image = qc_hdm.name;
%  im_show = imread(fullfile(subject.QCdir, image));
%  f = figure;
%  imshow(im_show);
%  uiwait(f)
 openfig(fullfile(subject.QCdir, image), 'visible')
 uiwait
 
 button=questdlg({'Headmodel okay?'},'Proceed with next subject','Yes');
 if strcmpi(button,'No') % discard subject if headmodel is bad
  dropout_file = fullfile(subject.analysis_dir, 'subjects_dropout.mat');
  load(dropout_file);
  reason = 'source problems';
  rej_subjects{end+1,1} = subject.id;
  rej_subjects{end,2} = reason;
  save('subjects_dropout.mat', 'rej_subjects')
 end
 
 close all
 clear subject
 
end

%% Select subjects for processing
% discard all subjects who have less than 30 trials and/or bad headmodel
% and move full subject-directory to a drop-out folder 

dropout_folder = '/home/uni10/nmri/projects/cstier/aging_processing/dropout';

for s = 1:length(rej_subjects(:,1))
 if exist(char(rej_subjects(s,1)), 'file')
 movefile(char(rej_subjects(s,1)), dropout_folder)
 end
end

proc_subjects = good_subjects;

% now do source projection and compute power and connectivity at sources

for i = 1:length(proc_subjects)
 
 % load subject struct
 subject = proc_subjects(i);
 subj_id = char(proc_subjects(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 load(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.mat'))
 subject = nmri_load_subject_most_advanced(subject);
 
 % Source analysis
 subject = cs_nmri_processing(subject);
 
end

%% Quality check of source results

for i = 1:length(proc_subjects)
 
 % load subject struct
 subject = proc_subjects(i);
 subj_id = char(proc_subjects(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 load(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.mat'))
 subject = nmri_load_subject_most_advanced(subject);

 button=questdlg({'Do you want to review source power plots now?'},'Review QC plot?','No');
 if strcmpi(button,'Yes')
  qc_power = dir(fullfile(subject.QCdir, 'pow_*'));
  for p = 1:length(qc_power)
   image = qc_power(p).name;
   im_show = imread(fullfile(subject.QCdir, image));
   f = figure;
   imshow(im_show);
   uiwait(f)
  end
 end

 button=questdlg({'Do you want to review source connectivity plots now?'},'Review QC plot?','No');
 if strcmpi(button,'Yes')
  qc_con = dir(fullfile(subject.QCdir, 'con_coh_img*'));
  for p = 1:length(qc_con)
   image = qc_con(p).name;
   im_show = imread(fullfile(subject.QCdir, image));
   f = figure;
   imshow(im_show);
   uiwait(f)
  end
 end
 
 button=questdlg({'source plots okay?'},'Proceed with next subject','Yes');
 if strcmpi(button,'No') % discard subject if headmodel is bad
  dropout_file = fullfile(subject.analysis_dir, 'subjects_dropout.mat');
  load(dropout_file);
  reason = 'sources results bad';
  rej_subjects{end+1,1} = subject.id;
  rej_subjects{end,2} = reason;
  save('subjects_dropout.mat', 'rej_subjects')
 end
 
 close all
 clear subject
 
end

%% filter again only clean and good subjects and save
if exist(dropout_file)
 load(dropout_file);
else
 rej_subjects = {};
 save('subjects_dropout.mat', 'rej_subjects')
 load(dropout_file);
end

if ~isempty(rej_subjects)
 all_subjects = setdiff(selected_subjects, rej_subjects(:,1)); 
 dropout_folder = '/home/uni10/nmri/projects/cstier/aging_processing/dropout';

 for s = 1:length(rej_subjects(:,1))
  if exist(char(rej_subjects(s,1)), 'file')
   if exist(dropout_folder)
    movefile(char(rej_subjects(s,1)), dropout_folder)
   else 
    mkdir dropout
    movefile(char(rej_subjects(s,1)), dropout_folder)
   end
  end
 end

else
 subjlist_final = selected_subjects;
end

% save list in mat-file
save('all_subjects_list.mat', 'subjlist_final');

% load subject-information and save for export of metrics and statistical analyses
% get subject arrays
all_subjects = {};

for i = 1:length(subjlist_final)
 
 % load subject struct
 subject = subjlist_final(i);
 subj_id = char(subjlist_final(i));
 exam_id = fullfile('ses-rest/meg', [subj_id '_ses-rest_task-rest_proc-sss.fif']);
 analysis_dir = pwd;
 load(fullfile(analysis_dir,subj_id, 'conf', 'ses-rest_meg_ses-rest_task-rest_proc-sss', 'subject_info.mat'))
 subject = nmri_load_subject_most_advanced(subject);
 
all_subjects{i,1} = subject;
 
end

save('all_subjects.mat', 'all_subjects'); 
% save('cohortX_final.mat', 'all_subjects'); % change file name accordingly in case if age subgroups will be analysed

 
   
   
   