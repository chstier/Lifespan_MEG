%% This scripts loads cortical thickness values for each subject and resamples them for functional resting-state networks
% (Yeo et al., 2011; https://doi.org/10.1152/jn.00338.2011)
%
% Written by Christina Stier, 2022

% load surface description for remapping from suma vertex-resolution to yeo-networks 
load('/home/uni10/nmri/projects/cstier/aging_analysis/conf/atlas/Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat','suma_all')

% load individual cortical thickness levels
metric = {'thickness_smoothed_12_'};
rootdir_th = ('/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/thickness/');

group = 'all';
yeo_all = {};

for var_metric =  metric{1}
       data = load([rootdir_th...
       metric{1}...
       'suma_all_N350.mat']); % change filenames here if needed
   for i = 1:length(data.thickness)
       m3 = data.thickness{i};
       m3(2005:2338,1) = NaN; % set subcortical nuclei to NaN
       regvalues_th = cs_nf_ft_report_regional(suma_all,m3);
       regvalues_th(1,:) = []; 
       roi = regvalues_th(:,1);
       regvalues_th(:,1:3) = [];
       yeo_all{i,1} = cell2mat(regvalues_th);
   end
end

   
opt.fwhm = '12';
filebase=fullfile(rootdir_th,[ 'thickness_smoothed_' opt.fwhm '_yeo_N' num2str(length(yeo_all))]);

all_subjects = data.all_subjects;
save([filebase '.mat'],'yeo_all','all_subjects')

nmri_write_mgh([filebase '.mgh'],eye(4),yeo_all)


%% make separate files for each network to run in PALM

load('/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/thickness/thickness_smoothed_12_yeo_N350.mat')
load('regnames_yeo_manual.mat')

mkdir /home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/thickness/yeo
report_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/thickness/yeo';

all = cell2mat(yeo_all');

hemi = {};
region = {};

for r = 1:length(roi)
 hemi{r,1} = roi{r,2}
 nw{r,1} = roi{r,3}
end
 

% load covariables
descr = csvread('/home/uni10/nmri/projects/cstier/aging_pipeline/design_example.csv'); % modify here. Should load co-variables for each subject

% concatenate cortical thickness values plus co-variables in a new
% design-file for statistical analysis
for i = 1:16
 network = (all(i,:))';
 network = [network descr]; 
 csvwrite(fullfile(report_dir, ['yeo7_N350_th_age_age2_sex_icv_', hemi{i}, '_', nw{i}, '.csv']), network)

end
 


 
 
 
 