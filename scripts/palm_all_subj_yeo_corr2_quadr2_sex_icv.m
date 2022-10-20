% Do PALM on IGE - controls (ImCoh/Power ~ group)

analysis_dir = '/home/uni10/nmri/projects/cstier/aging_analysis/export/all_subjects';
export_func = '/home/uni10/nmri/projects/cstier/aging_analysis/export/all_subjects/yeo/';
export_th = '/home/uni10/nmri/projects/cstier/aging_analysis/export/all_subjects/ld10_sm/yeo/';
results_dir = '/home/uni10/nmri/projects/cstier/aging_analysis/results/yeo_func_struct_age_age2_sex_icv';

% palm -i export/all_subjects/yeo/7networks_yeo_N350_coh_img_Alpha_left_default.csv -d export/all_subjects/ld10_sm/yeo/yeo7_N350_th_age_age2_sex_icv_left_default.csv -demean -t con_yeo_th_age_age2_icv_sex.csv -saveglm -savedof -accel tail -n 500 -o results/yeo_func_struct_age_age2_sex_icv/Alpha/coh_img


load('regnames_yeo_manual.mat')
hemi = {};
region = {};

for r = 1:length(roi)
 hemi{r,1} = roi{r,2}
 nw{r,1} = roi{r,3}
end

%fullfile(report_dir, ['yeo7_N350_th_age_age2_sex_icv_', hemi{i}, '_', nw{i}, '.csv']), network)
 

% loop over frequencies
% loop over metric
% loop over networks

metric = {'coh_img_', 'power_'};
freq = {'Delta','Theta','Alpha','Beta1','Beta2','Gamma'};
scale = {'abs_not_scaled_yeo_'};

for var_freq = 1:length(freq)
 for var_metric = 1:length(metric)
  for i = 1:16
   func_file = fullfile(export_func, ['7networks_yeo_N350_', metric{var_metric}, freq{var_freq}, '_', hemi{i}, '_', nw{i}, '.csv']);
   design_file = fullfile(export_th, ['yeo7_N350_th_age_age2_sex_icv_', hemi{i}, '_', nw{i}, '.csv']);
  
   results = fullfile(results_dir, freq{var_freq}, [metric{var_metric}, hemi{i}, '_', nw{i}]);
  
   code_line = [{'-i', func_file, '-d', design_file, '-demean', '-t', 'con_yeo_th_age_age2_icv_sex.csv', '-saveglm', '-savedof', '-accel', 'tail', '-n', '500', '-o', results}];
  
   palm(code_line{:})
  end
 end
end

%   palm -i func_file -d design_file -demean -t con_yeo_th_age_age2_icv_sex.csv -saveglm -savedof -accel tail -n 500 -o results


% %coh_img, power, wpli_debiased for Alpha
palm -i reports/all_subjects/global_abs_mean_coh_img_Alpha_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Alpha/coh_img
palm -i reports/all_subjects/global_abs_mean_power_Alpha_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Alpha/power

%coh_img, power, wpli_debiased for Beta1
palm -i reports/all_subjects/global_abs_mean_coh_img_Beta1_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Beta1/coh_img
palm -i reports/all_subjects/global_abs_mean_power_Beta1_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Beta1/power

%coh_img, power, wpli_debiased for Beta2
palm -i reports/all_subjects/global_abs_mean_coh_img_Beta2_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Beta2/coh_img
palm -i reports/all_subjects/global_abs_mean_power_Beta2_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Beta2/power

%coh_img, power, wpli_debiased for Delta
palm -i reports/all_subjects/global_abs_mean_coh_img_Delta_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Delta/coh_img
palm -i reports/all_subjects/global_abs_mean_power_Delta_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Delta/power

%coh_img, power, wpli_debiased for Gamma
palm -i reports/all_subjects/global_abs_mean_coh_img_Gamma_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Gamma/coh_img
palm -i reports/all_subjects/global_abs_mean_power_Gamma_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Gamma/power

%coh_img, power, wpli_debiased for Theta
palm -i reports/all_subjects/global_abs_mean_coh_img_Theta_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Theta/coh_img
palm -i reports/all_subjects/global_abs_mean_power_Theta_values.csv -d design_all_corr2_sex_manualdemean_icv.csv -demean -t con_corr2_quadr2_all_sex_icv.csv -saveglm -savedof -accel tail -n 500  -o results/global_corr2_quadr2_sex_icv_manualdemean_plusdemean/Theta/power


% get overview results

freqname = {'Delta_', 'Theta_', 'Alpha_', 'Beta1_', 'Beta2_', 'Gamma_'};
freq = {'Delta', 'Theta', 'Alpha', 'Beta1', 'Beta2', 'Gamma'};
metric = {'coh_img_', 'power_'};
analysis_type = {'dat_'};
analysis_type_2={'tstat_', 'tstat_fwep_'};
contrasts = {'c1', 'c2', 'c3', 'c4', 'c5', 'c6'};

stat_table = {};
stat_head = {'frequency', 'metric', 'contrasts', 't', 'p(fwe)'};
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
%                           case 'p(uncp)'
%                               stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{2} contrasts{1, var_contrasts} '.csv']))); 
                          case 'p(fwe)'
                              stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{2} contrasts{1, var_contrasts} '.csv'])));   
%                           case 'p(fwe)_eeg'
%                               stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{4} contrasts{1, var_contrasts} '.csv'])));    
%                           case 't_meg'
%                               stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{5} contrasts{1, var_contrasts} '.csv'])));                      
%                           case 'p(fwe)_meg'
%                               stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{6} contrasts{1, var_contrasts} '.csv'])));    
%                           case 'cohens_d_eeg'
%                               stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{7} contrasts{1, var_contrasts} '.csv']))); 
%                           case 'cohens_d_meg'
%                               stat_table{stat_row,pl} = csvread(char(fullfile(rootdir, [metric{1, var_metric} analysis_type{1} analysis_type_2{8} contrasts{1, var_contrasts} '.csv']))); 
                      end
                  end
              stat_row=stat_row+1;         
        end
    end
end

nf_csvwrite(fullfile(results_dir,['global_abs_mean_results_' 'all_subj_corr2_quadr2_sex_icv_500_tail_manualdemean_plusdemean_pvalue' '.csv']),stat_table);


