%% This scripts runs the spin-test for the correlation between the age-effects on cortical thickness and MEG levels
% Based on the toolbox by Alexander-Bloch et al., 2018: https://doi.org/10.1016/j.neuroimage.2018.05.070
%
% Step 1: 
% Load the statistical map for the age effects on cortical thickness (here t-stats for decrease (C2) or quadratic effect (C4))
% Spin data for the left and right hemisphere seperately (cs_SpinPermuFS)

% Step 2:
% Load the statistical maps for the age effeccts on MEG markers and compute the real Spearman's correlation with the statistical map 
% for cortical thickness. Test Spearman's correlation coefficient against the null distribution created using SpinPermuFS and save output
% 
% written and edited by Christina Stier, 2022

%% STEP 1 
% prepare all the files and folders needed
suma_map = fullfile(pwd, 'conf/suma-all-fsaverage-10.mat'); % get SUMA incl. subcortical nuclei (ld 10)
suma_all_all = load(suma_map);

suma_cortex = fullfile(pwd, 'conf/suma-cortex-fsaverage-10.mat'); % get SUMA without subcortical nuclei
load(suma_cortex);

mkdir 'Rotation_corrmaps'

% get thickness contrast (here t-statistics dpv as output from PALM-analysis) and assign to left and right hemisphere separately 
stats_dir_th = fullfile(pwd, 'results/thickness_age_age2_sex_tiv'); % change name of results-folder if needed
mod2 = 'dpv_'; % choose modality (here: results per-vertex)
th_file = '_dpv_tstat_c4.mgz'; % change to contrast needed here (linear / quadratic effects)

th_t = load_mgh(fullfile(stats_dir_th, th_file));
th_t = th_t(1:2004,:); % only use cortical nuclei

% apply mask (to make sure that all medial points are NaN) 
th_t(suma_all.msk == 0) = NaN;

th_left = th_t(suma_all.hemi == 1,:); % subcortical values are already zeroed out
th_right = th_t(suma_all.hemi == 2,:);

readleft = th_left;
readright = th_right;

% choose number of permutations:
permno = 10000;
perm = '10000';

% spin and generate 0 distribution

wsname = fullfile(pwd, 'Rotation_corrmaps', ['thickness_nosubc_dpv_tstat_c4_' perm 'perm_spearman.mat']); % change filename of output if needed

cs_SpinPermuFS(readleft, readright, permno, wsname) 

%% STEP 2
% indicate left and right for the second modality
freqname = {'Delta','Theta','Alpha', 'Beta1', 'Beta2', 'Gamma'};
metric_name = {'coh_img', 'power'};
mod = {'_dpv_tstat_c1', '_dpv_tstat_c2', '_dpv_tstat_c3', '_dpv_tstat_c4'};

stats_dir_func = fullfile(pwd, 'results/MEG_age_age2_sex_tiv'); % change name of results-folder if needed 
results_dir = fullfile(pwd, 'Rotation_corrmaps');

for f = 1:length(freqname)
 for m = 1:length(metric_name) 
  for o = 1:length(mod)
   
   func_file = fullfile(stats_dir_func, freqname{f}, [metric_name{m} mod{o} '.mgz']);
   func_name = [freqname{f} '_' metric_name{m} mod{o}];
   func_t = load_mgh(func_file);
   
   func_t = func_t(1:2004,:); % delete values for subcortical nuclei
   
   % apply mask (to make sure that all medial points are NaN)
   func_t(suma_all.msk == 0) = NaN;

   [pval, realrho, realpval, nullroh] = cs_pvalvsNull_spearman(th_t, func_t, permno,wsname); 
   
   % save results in struct and csv
   stats(o).pval = pval;
   stats(o).real_rho = realrho;
   stats(o).realpval = realpval;
   stats(o).mod = func_name;
   
   filename_stats = fullfile(results_dir, ['corr_thickn_c4_' freqname{f} '_' metric_name{m} '_' perm '_perm_spearman.csv']);
   
   writetable(struct2table(stats), filename_stats);
 
  end
 end
end

%% concat all results into one file for overview
rootdir = results_dir;

for m = 1:length(metric_name)  
  
  out1 = readtable(fullfile(results_dir, ['corr_thickn_c4_' freqname{1} '_' metric_name{m} '_' perm '_perm_spearman.csv'])); % change contrast if needed
  out2 = readtable(fullfile(results_dir, ['corr_thickn_c4_' freqname{2} '_' metric_name{m} '_' perm '_perm_spearman.csv']));
  out3 = readtable(fullfile(results_dir, ['corr_thickn_c4_' freqname{3} '_' metric_name{m} '_' perm '_perm_spearman.csv']));
  out4 = readtable(fullfile(results_dir, ['corr_thickn_c4_' freqname{4} '_' metric_name{m} '_' perm '_perm_spearman.csv']));
  out5 = readtable(fullfile(results_dir, ['corr_thickn_c4_' freqname{5} '_' metric_name{m} '_' perm '_perm_spearman.csv']));
  out6 = readtable(fullfile(results_dir, ['corr_thickn_c4_' freqname{6} '_' metric_name{m} '_' perm '_perm_spearman.csv']));
  
  allres = [out1; out2; out3; out4; out5; out6];
  
  filename_allres = fullfile(results_dir, ['corr_thickn_c4_' mod2 metric_name{m} '_allfreqs_' perm '_perm_spearman.csv']);
   
  writetable(allres, filename_allres);
  
  clear allres
  
end
