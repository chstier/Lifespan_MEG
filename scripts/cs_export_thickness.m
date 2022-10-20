%% This script creates cortical thickness maps for each subject and exports for statistical analyses
% Written and modified by Christina Stier, 2020/22

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

% orelse load all_subjects from previous analyses
analysis_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline';
outputdir = '/export/cohort_paper/thickness';
load('all_subjects.mat') % orelse if only specific cohort is needed: load('cohortX_final.mat');

%% export vertex-wise
thickness = {};

% get thickness-files
opt = [];

% smooth before
opt = [];
opt.sigma = '2'; % choose sigma
opt.fwhm = '12'; % choose FWHM
opt.ld = '10'; % choose SUMA resolution (ld-factor)
opt.subc = '334'; % indicate how many subcortical nuclei we have (depending on the resolution). Will be set to 0.

suma_folder = '/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/'; % change accordingly
mod = '_anat_T2_T1'; % change to your needs

for i = 1:length(all_subjects)
 subject = all_subjects{i,1};
 
[smoothed_right, smoothed_left] = cs_smooth_thickness_suma_all(subject, opt, suma_folder, mod);
 
 opt.right_file = smoothed_right;
 opt.left_file = smoothed_left;
 cu_suma_all = cs_get_thickness_suma_all_sm(subject, opt, suma_folder, mod);
 thickness{i,1} = cu_suma_all;

end

%% save output in .mat file and write mgh-file for palm-analyses with positive

% regular values
mkdir([analysis_dir, outputdir]);
filebase=fullfile(analysis_dir, outputdir,[ 'thickness_smoothed_' opt.fwhm '_suma_all_N' num2str(length(all_subjects))]);

save([filebase '.mat'],'thickness','all_subjects')
nmri_write_mgh([filebase '.mgh'],eye(4),thickness)

%% export global metrics (saves .csv-file in reports-folder)
opt = [];
opt.metrics = 'thickness';
opt.statistics = 'mean';
opt.sigma = '2';
opt.fwhm = '12';
opt.ld = '10';
opt.subc = '334'; % 334 in case of ld 10; 5390 in case of ld 40
ld = '10';

stat_items = {'Subject_ID','Exam_ID','Metric','Stat','Global'};
stat_cols=length(stat_items); 

stat_out=cell((length(all_subjects))+1,stat_cols);

% Print Header
stat_out(1,1:length(stat_items))=stat_items;
stat_row=2;

% load export file in which vertex-wise values are stored
filebase=fullfile(analysis_dir, outputdir,[ 'thickness_smoothed_' opt.fwhm '_suma_all_N' num2str(length(all_subjects))]);
load(filebase)

for i = 1:length(thickness)
 subject = thickness{i,1};
 
 % set subcortical nuclei to NaN and compute global value
 subject(2005:2338) = NaN;
 subject_global = nanmean(subject);
%  global_th{i,1} = subject_global; 
 
 % put default stat_items
 for si=1:length(stat_items)
  switch stat_items{si}
   case 'Subject_ID'
    stat_out{stat_row,si}=all_subjects{i}.id;
    case 'Exam_ID'
    stat_out{stat_row,si}=all_subjects{i}.exam_id;
   case 'Metric'
    stat_out{stat_row,si}=opt.metrics;
 %   case 'Freq'
 %    stat_out{stat_row,si}=freq_txt{f};       
   case 'Stat'
    stat_out{stat_row,si}=opt.statistics;
   case 'Global'
    stat_out{stat_row,si}=subject_global;
  end 
 end
 % next subject
 stat_row=stat_row+1;
end
 
% we always report
opt.report=true;

if ~isfield(opt,'report') || ~islogical(opt.report)
 opt.report=false;
else
 % we want reports
 if ~isfield(opt,'report_dir')
  opt.report_dir=fullfile(pwd,'reports/cohort_paper/thickness/'); % change folder if needed
 end
 if ~isfield(opt,'report_regional') || ~islogical(opt.report_regional)
  opt.report_regional=false;
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
end

% only save global values for global group analyses
stat_value = stat_out(2:end,end);
% save as csv
nf_csvwrite(fullfile(opt.report_dir,['global_' opt.statistics '_' opt.metrics '_smoothed_' opt.fwhm '_ld' ld '_descr.csv']),stat_out); % save full descr
nf_csvwrite(fullfile(opt.report_dir,['global_' opt.statistics '_' opt.metrics '_smoothed_' opt.fwhm '_ld' ld '_values.csv']),stat_value); % sa



