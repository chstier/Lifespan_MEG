%% This scripts loads MEG connectivity and power values for each subject and resamples them for functional resting-state networks
% (Yeo et al., 2011; https://doi.org/10.1152/jn.00338.2011)
%
% Written by Christina Stier, 2022

% load surface description for remapping from suma vertex-resolution to yeo-networks 
load('/home/uni10/nmri/projects/cstier/aging_analysis/conf/atlas/Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat','suma_all')

% load individual MEG levels
metric = {'coh_img_', 'power_'};
freq = {'Delta_','Theta_','Alpha_','Beta1_','Beta2_','Gamma_'};
scale = {'abs_not_scaled_'};
rootdir = ('/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/MEG/');

group = 'all';
yeo_all = {};

for var_freq = 1:length(freq)
 for var_metric = 1:length(metric)
            data = load([rootdir...
            metric{var_metric}...
            freq{var_freq}...
            scale{1}...
            'N350.mat']); % change filenames here if needed
        for i = 1:length(data.scale_metrics) 
            m3 = data.scale_metrics{i};
            m3(2005:2338,1) = NaN; % set subcortical nuclei to NaN
            regvalues_th = cs_nf_ft_report_regional(suma_all,m3);
            regvalues_th(1,:) = []; 
            roi = regvalues_th(:,1);
            regvalues_th(:,1:3) = [];
            yeo_all{i,1} = cell2mat(regvalues_th);
        end
        
   filebase=fullfile(rootdir, 'yeo', [ metric{var_metric} freq{var_freq} scale{1} 'yeo_N' num2str(length(yeo_all))]);
   all_subjects = data.all_subjects;
   save([filebase '.mat'],'yeo_all','all_subjects')

   nmri_write_mgh([filebase '.mgh'],eye(4),yeo_all)     
 end
end

 
%% make separate files for each network to run in PALM
load('regnames_yeo_manual.mat')

hemi = {};
region = {};

for r = 1:length(roi)
 hemi{r,1} = roi{r,2}
 nw{r,1} = roi{r,3}
end
 
report_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/MEG/yeo';
mkdir /home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/MEG/yeo
rootdir = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/yeo/';

metric = {'coh_img_', 'power_'};
freq = {'Delta_','Theta_','Alpha_','Beta1_','Beta2_','Gamma_'};
scale = {'abs_not_scaled_yeo_'};

for var_freq = 1:length(freq)
 for var_metric = 1:length(metric)
            data = load([rootdir...
            metric{var_metric}...
            freq{var_freq}...
            scale{1}...
            'N350.mat']); % change filenames here if needed
           
     all = cell2mat(data.yeo_all');
     
   for i = 1:16
    network = (all(i,:))';
    csvwrite(fullfile(report_dir, ['7networks_yeo_N350_', metric{var_metric}, freq{var_freq}, hemi{i}, '_', nw{i}, '.csv']), network)
   end
  
 end
end
