%% This script exports selected variables (power/connectivity) for further statistical comparisons between groups
% It will create a folder (export), in which .mat and .mgh files with vertex-values will be
% stored, and a .csv-file in the /reports-folder containing global values
% for each frequency band and metric

% Written by Christina Stier 2020/22

%% If not done yet, load final subject-list and load subject-struct for each subject and save for export of metrics
load('all_subjects_list.mat')

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

%% Export metrics for all subjects as used in the current manuscript (Stier et al.)
% needed for statistical analysis of vertex-based connectivity and power

load('all_subjects.mat')
cohort_stat = 'all_subjects';
analysis_dir = pwd;

% export metrics to mgh-file for each individual
opt=[];
opt.metrics = {'coh_img','power'}; % change to your needs
opt.scale={'none'};
opt.global={'abs'};
opt.dim_reduction={'mean'};
opt.save_average={'true'};
opt.output = [analysis_dir, '/export/cohort_paper/MEG/'];
nmri_export_metrics(all_subjects,opt)

% generate tabulated overview for global metrics (averaged across the individual brain)
opt.statistics={'abs_mean'};
opt.report_dir = [analysis_dir, '/reports/cohort_paper/MEG/'];
[stat_out] = cs_create_global_files(all_subjects, opt);

%% Export metrics for respective age groups
% needed for visualization of global power and connectivity in R

analysis_dir = pwd;
cohort = {'cohortX'}; % {'cohort1', 'cohort2','cohort3', 'cohort4', 'cohort5', 'cohort6', 'cohort7'};

for c = 1:length(cohort)
 
  % load cohort-files containing subject-infos
  cohort_file = [cohort{c} '_final.mat'];
  load(cohort_file);
  
  cohort_stat = cohort(c);
  
  if ~exist('/export')
   mkdir export
  end
  
  % export metrics to mgh-file for each individual and cohort
  opt=[];
  opt.metrics = {'coh_img','power'}; % change to your needs
  opt.scale={'none'};
  opt.global={'abs'};
  opt.dim_reduction={'mean'};
  opt.save_average={'true'};
  opt.output = fullfile(analysis_dir, 'export', cohort{c});
  nmri_export_metrics(all_subjects,opt)
  
  % generate tabulated overview for global metrics (averaged across the individual brain)
  opt.statistics={'abs_mean'};
  opt.report_dir = [analysis_dir, '/reports/cohort_paper/MEG/'];
  [stat_out] = cs_nmri_statistics_aging(all_subjects,cohort_stat, opt);
  
  clear all_subjects cohort_stat

end
  
