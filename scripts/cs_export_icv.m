% This script loads subject-list of interest, extract intracranial volumes for each
% subject and saves in a csv-file
% Christina Stier, written and modified in 2020/22

outputdir = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/icv';

% create folder, in which head size values should be stored
if ~exist(outputdir, 'dir')
 mkdir ./export/cohort_paper icv
end

suma_folder = '/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/';
mod = '_anat_T2_T1';

load('all_subjects.mat') % orelse load('cohortX_final.mat') if only separate cohorts are tested
subjlist = 'all_subjects';

icv_all = {};

for i = 1:length(all_subjects)
 subject = all_subjects{i,1};
 icv = cs_get_icv(subject, suma_folder, mod);
 
 icv_all{i,1} = icv;
 clear icv
end

save(fullfile(outputdir,['icv_' subjlist '.mat']),'icv_all');
dlmwrite(fullfile(outputdir,['icv_' subjlist '.csv']), icv_all, 'delimiter', ',', 'precision', 8);
 
