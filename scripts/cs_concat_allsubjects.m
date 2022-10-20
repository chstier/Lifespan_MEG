%% This script concatenates the subject infos of several cohorts

analysis_dir = pwd;

cohort = {'cohort1', 'cohort2','cohort3', 'cohort4', 'cohort5', 'cohort6', 'cohort7'};


all_cohorts = {};

for c = 1:length(cohort)

  % load cohort-files containing subject-infos
 
  cohort_file = [cohort{c} '_final.mat'];
  load(cohort_file);
    
  s{c}.cohort = all_subjects;  

end

all_subjects = [s{1}.cohort; s{2}.cohort; s{3}.cohort; s{4}.cohort; s{5}.cohort; s{6}.cohort; s{7}.cohort]

save('all_subjects.mat', 'all_subjects') 

% all_subjects = [s{1}.cohort; s{2}.cohort; s{3}.cohort; s{4}.cohort; s{5}.cohort]
% all_subjects = [s{1}.cohort; s{2}.cohort]
  
%   % export metrics to mgh-file for each individual and cohort
%   opt=[];
%   opt.metrics = {'coh_img','power'}; % already present: 'power' {'coh_img', 'wpli_debiased'}
%   opt.scale={'none'};
%   opt.global={'abs'};
%   opt.dim_reduction={'mean'};
%   opt.save_average={'true'};
%   opt.output= fullfile(analysis_dir, 'export', 'cohorts6_7')
%   nmri_export_metrics(all_subjects,opt)
% 
%  
%   % export for graph analysis
%   
%   opt=[];
%   opt.metrics = {'coh_img'}; % already present: 'power' {'coh_img', 'wpli_debiased'}
%   opt.scale={'none'};
%   opt.global={'abs'};
%   opt.dim_reduction={'none'};
%   opt.save_mgh=false; % no mgh
%   opt.output= fullfile(analysis_dir, 'GTA_export', 'all_subjects');
%   nmri_export_metrics(all_subjects,opt)
% 
%   
% 
% 
% 
