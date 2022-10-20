% This script describes the statistical analysis of age effects on MEG power and connectivity as well as cortical thickness. 
% Steps for visulatization of the results are described and are partially done using Matlab or R.
% 
% Christina Stier 2022

%% Vertex-wise analysis using PALM (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/UserGuide)
% Reference: Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation inference for the general linear model. NeuroImage, 2014;92:381-397 

% Preparations:
% set up regressors: get intracranial volume if needed using cs_export_icv.m;
% prepare design and contrast files for the statistical analysis (PALM) and save: see example for the sample used in the current manuscript (Stier et al.)

% example file for the statistical contrasts: contrasts_example.csv: 
% tests the positive and negative linear effects of age (C1 and C2), the quadratic (U- and inverted U-shaped) effects of age (C3 and C4), the effects of sex (C5) and intracranial volume (C6). 
% example file for the design as used in the current manuscript (Stier et al.): design_example.csv

% example analysis for age-effects on MEG connectivity
analysis_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline';
results_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/results/MEG_age_age2_sex_tiv';

palm -i export/cohort_paper/MEG/coh_img_Beta1_abs_not_scaled_N350_all.mgh -m export/cohort_paper/all_suma_msk.mgh -s conf/suma-all-fsaverage-10.gii -d design_example.csv -demean -t contrasts_example.csv -T -tfce2d -logp -saveglm -savedof -accel tail -n 500 -o results/MEG_age_age2_sex_tiv/Beta1/coh_img

% example for visualization of the results (MEG; p-value):
% single plots for PALM results on low beta (contrasts 1-6)
load([analysis_dir '/conf/suma-all-fsaverage-10.mat'],'suma_all')

opt=[];
opt.per_hemi=1;
opt.per_cortex=1;
opt.rot=[90 0 ; -90 0];
opt.thresh=1.3;
opt.clim=[1.3 3.8]; % set thresholds accordingly (will show -log10(p))
opt.colormap='hot';
opt.colorbar='hot';
opt.scale=1;

freqname = {'Beta1_'};
metric = {'coh_img_'}; % change to the respective variables {'power_'}
analysis_type = {'tfce_'};
analysis_type_2={'tstat_fwep_', 'tstat_uncp_'};
contrasts = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6'}; 

rootdir = [results_dir '/Beta1/'];

for var_metric = 1:length(metric)
    for var_analysis_type = 1:length(analysis_type)
        for var_analysis_type_2 = 1:length(analysis_type_2)
            for var_contrasts = 1:length(contrasts)
                
                log_p_map=load_mgh([rootdir...
                    metric{1, var_metric}...
                    analysis_type{1,var_analysis_type} analysis_type_2{1, var_analysis_type_2}...
                    contrasts{1, var_contrasts} '.mgz']);

                opt.title=strcat(cell2mat(freqname),cell2mat(metric(1, var_metric)),cell2mat(analysis_type(1,var_analysis_type)),cell2mat(analysis_type_2(1, var_analysis_type_2)),cell2mat(contrasts(1, var_contrasts))) ;
                opt.output=strcat(fig_dir,cell2mat(metric(1, var_metric)),cell2mat(analysis_type(1,var_analysis_type)),cell2mat(analysis_type_2(1, var_analysis_type_2)),cell2mat(contrasts(1, var_contrasts)),'.png');
                
                hFig = nmri_plot_surface_suma(suma_all, log_p_map, opt);
            end
        end
    end
end

% example analysis for age-effects on cortical thickness
palm -i export/cohort_paper/thickness/thickness_smoothed_12_suma_all_N350.mgh -m export/cohort_paper/all_suma_msk_zerosubc.mgh -s conf/suma-all-fsaverage-10.gii -d design_example.csv -demean -t contrasts_example.csv -T -tfce2d -logp -saveglm -savedof -accel tail -n 500 -o results/thickness_age_age2_sex_tiv/

% use the example code above for the visualization of the results for cortical thickness
% script to make plots as in Stier et al. described (t-statistics for cortical thickness: decrease and inverted U-shape with age for significant regions:)
% cs_plot_sign_effects_palm.m

%% Global analysis using PALM (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/UserGuide)
% Reference: Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation inference for the general linear model. NeuroImage, 2014;92:381-397 

% example analysis for age-effects on global MEG levels (averaged markers across the individual brain)
palm -i reports/cohort_paper/MEG/global_abs_mean_coh_img_Beta1_values.csv -d design_example.csv -demean -t contrasts_example.csv -saveglm -savedof -accel tail -n 500  -o results/global_MEG_age_age2_sex_tiv/Beta1/coh_img

% make tables for results
results_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/results/global_MEG_age_age2_sex_tiv';

freqname = {'Beta1_'}; % modify {'Delta_', 'Theta_', 'Alpha_', 'Beta1_', 'Beta2_', 'Gamma_'}
freq = {'Beta1'}; % modify {'Delta', 'Theta', 'Alpha', 'Beta1', 'Beta2', 'Gamma'}

metric = {'coh_img_'}; % modify {'coh_img_', 'power_'}
analysis_type = {'dat_'};
analysis_type_2={'tstat_', 'tstat_fwep_'}; 
contrasts = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6'};

stat_table = {};
stat_head = {'frequency', 'metric', 'contrasts', 't', 'p'};
stat_table(1,1:length(stat_head)) = stat_head;
stat_row = 2;

for i=1:length(freqname)
    rootdir = fullfile(results_dir, '/', freq{i}, '/');
    for var_metric = 1:length(metric)
            for var_contrasts = 1:length(contrasts)

                  for pl = 1:length(stat_head)
                      switch stat_head{pl}
                          case 'frequency'
                              stat_table{stat_row,pl} = freq{i};
                          case 'metric'
                              stat_table{stat_row,pl} = metric{1, var_metric};
                          case 'contrasts'
                              stat_table{stat_row,pl} = contrasts{1, var_contrasts};
                          case 't'
                              stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{1} contrasts{1, var_contrasts} '.csv'])));  
                          case 'p'
                              stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{2} contrasts{1, var_contrasts} '.csv'])));   
                      end
                  end
              stat_row=stat_row+1;         
        end
    end
end

nf_csvwrite(fullfile(results_dir,['global_abs_mean_results_' 'all_subj_age_age2_sex_tiv' '.csv']),stat_table);

% to compute partial correlation coefficients based on the t-statistics, use the script 'cs_partialr_palm.m'

% for plotting the data distributions of global values in R:  
% get the following files and save in your R analysis folder:
% - 'cohortX_metric_summary..'-file (global MEG values) from /reports/MEG-folder
% - 'metric_summary_...'-file from /reports/thickness-folder containing the global cortical thickness values
% - 'Demographics_cohortX_final.mat'-file from /reports-folder
% - 'icv_all_subjects.csv'-file from /exports/cohortX/icv-folder and
% create datasets in R using the scripts 'global_datasets_...R'
% plot distributions with age using the scripts 'global_plots...R'

%% Correlation analysis hetween the age effects on cortical thickness and MEG markers using the spin-test 

% find the scripts for these steps in /scripts/utilities/cs_spatialpermutation_test
% run cs_spin_corrmaps.m and find results in /pwd/Rotation_corrmaps
% the analysis is based on the spin-test toolbox by Alexander-Bloch et al., 2018: https://doi.org/10.1016/j.neuroimage.2018.05.07

% if you like to visualize the correlation in R, copy the results files into your R analysis folder
% (age effects (t-values) on MEG markers for each frequency band: e.g. Delta_coh_img_dpv_tstat_c1.mat and those for cortical thickness e.g. _dpv_tstat_c2.mat)
% then use the R script 'corr_t_struct_func.R'

%% Regression analysis for cortical thickness on MEG levels across the adult lifespan sample

% the analysis was carried out for 7 functional resting-state networks (Yeo et al., 2011; https://doi.org/10.1152/jn.00338.2011).
% do resampling of MEG/cortical thickness values for the yeo-networks first:
% use cs_export_yeo_function.m and cs_export_yeo_thickness.m

% Run PALM analyis for each network. Prepare folders.
analysis_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper';
export_func = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/MEG/yeo';
export_th = '/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohort_paper/thickness/yeo';
results_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/results/yeo_func_struct_age_age2_sex_tiv';

% load network names
load('regnames_yeo_manual.mat')
hemi = {};
region = {};

for r = 1:length(roi)
 hemi{r,1} = roi{r,2};
 nw{r,1} = roi{r,3};
end

% create a contrast file, which tests positive (c1) and negative prediction (c2) of
% cortical thickness on MEG levels (while regressing out the effects of age, age2, sex, tiv)
con = [1 0 0 0 0 ; -1 0 0 0 0];
csvwrite('con_yeo_th_age_age2_tiv_sex.csv', con)

% Run PALM in a loop
metric = {'coh_img_', 'power_'};
freq = {'Delta','Theta','Alpha','Beta1','Beta2','Gamma'};
scale = {'abs_not_scaled_yeo_'};
n_subj = {'N350_'};

% loop over frequencies
for var_freq = 1:length(freq)
 % loop over metric
 for var_metric = 1:length(metric)
  % loop over networks
  for i = 1:16
   func_file = fullfile(export_func, ['7networks_yeo_', n_subj{1}, metric{var_metric}, freq{var_freq}, '_', hemi{i}, '_', nw{i}, '.csv']); % modify N here
   design_file = fullfile(export_th, ['yeo7_', n_subj{1}, 'th_age_age2_sex_icv_', hemi{i}, '_', nw{i}, '.csv']);
  
   results = fullfile(results_dir, freq{var_freq}, [metric{var_metric}, hemi{i}, '_', nw{i}]);
  
   code_line = [{'-i', func_file, '-d', design_file, '-demean', '-t', 'con_yeo_th_age_age2_tiv_sex.csv', '-saveglm', '-savedof', '-accel', 'tail', '-n', '500', '-o', results}];
  
   palm(code_line{:})
  end
 end
end

% concatenate the results for overview and plotting in R
hemi(1) = []; % exclude medial wall
hemi(8) = []; 
nw(1) = [];
nw(8) = [];

analysis_type = {'dat_tstat_', 'dat_tstat_fwep_'};
contrasts = {'c1', 'c2'}; % get results for positive/negative prediction by thickness

stat_table = {};
stat_head = {'frequency', 'network', 'contrasts', 't', 'p'};
stat_table(1,1:length(stat_head)) = stat_head;
stat_row = 2;

for m = 1:length(metric)
 for c = 1:length(contrasts) 
   for f = 1:length(freq)
    for ii = 1:14
     
     for pl = 1:length(stat_head)
      switch stat_head{pl}
        case 'frequency'
          stat_table{stat_row,pl} = freq{f};
        case 'network'
          stat_table{stat_row,pl} = [hemi{ii}, '_', nw{ii}];
        case 'contrasts'
          stat_table{stat_row,pl} = contrasts{1};
        case 't' 
          stat_table{stat_row,pl} =  csvread(fullfile(results_dir, freq{f}, [metric{m}, hemi{ii}, '_', nw{ii}, '_', analysis_type{1}, contrasts{c} '.csv']));
        case 'p'  
          stat_table{stat_row,pl} =  csvread(fullfile(results_dir, freq{f}, [metric{m}, hemi{ii}, '_', nw{ii}, '_', analysis_type{2}, contrasts{c} '.csv']));
       end
     end 
      
      stat_row=stat_row+1;
      
    end
   end
   nf_csvwrite(fullfile(results_dir,['yeo_func_struct_allfreq_' metric{m} contrasts{c} '.csv']),stat_table);
 end
end

% for visualization of the results for the yeo networks using R, copy
% output into your R analysis folder ('yeo_func_struct_allfreq_cohimg/power_c1/2.csv') 
% and use R scripts: 'visualize_palm_corr_func_struct_yeo_imcoh.R' and 'visualize_palm_corr_func_struct_yeo_power_reverse.R'  

%% Subanalysis for differing trajectories between males and females (interaction effects between age and sex) using PALM (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/UserGuide)
% Reference: Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation inference for the general linear model. NeuroImage, 2014;92:381-397 

% Preparations
% prepare design and contrast files for the statistical analysis (PALM) and save
design_example = [0	1	0	-29	0	841	0	32456; 1	0	-29	0	841	0	-40673	0]; % for 2 subjects: columns 1 and 2 correspond to sex, columns 3 and 4 to age (demeaned), columns 5 and 6 to (demeaned age)^2 , columns 7 and 8 to total intracranial volume 
contrast_example = [0	0	-1	1	0	0	0	0; 0	0	1	-1	0	0	0	0; 0	0	0	0	-1	1	0	0; 0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	-1	1; 0	0	0	0	0	0	1	-1]; % contrasts 1-4 test whether there are differences between the sexes in the slopes

% example statistical analysis for interaction effects of sex and age on MEG signals
analysis_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline';
results_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/results/interaction_sex_age_age2_tiv';

palm -i export/cohort_paper/MEG/coh_img_Beta1_abs_not_scaled_N350_all.mgh -m export/cohort_paper/all_suma_msk.mgh -s conf/suma-all-fsaverage-10.gii -d design_interaction_sex_age_age2_tiv.csv -demean -t con_interaction_sex_age_age2_tiv.csv -T -tfce2d -logp -saveglm -savedof -accel tail -n 500 -o results/interaction_sex_age_age2_tiv/Beta1/coh_img

% use example code from main vertex-analysis for plotting the results