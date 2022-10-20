% This scripts converts t-values into partial correlatin coefficients 
% (see e.g.  R. Rosenthal, H. Cooper, L. Hedges, Parametric measures of effect size. The handbook of research synthesis 621, 231-244 (1994))

% Formula
% r = sign(t)*sqrt(t^2 / N - rank(M) + t^2)

% Load the results of the models (MEG markers ~ age + age^2 + sex + tiv analysed using PALM)
% load t-values for each line of the table

M = [1 0 0 0; -1 0 0 0; 0 1 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]; % indicate design matrix
N = 350; % number of subjects
tbl = readtable('/home/uni10/nmri/projects/cstier/aging_pipeline/results/global_MEG_age_age2_sex_tiv/global_abs_mean_results_all_subj_age_age2_sex_tiv.csv');
t_values = tbl.t;


r_partial = {};
R2_partial = {};

for tv = 1:length(t_values)
 t = t_values(tv);
 r = sign(t)*sqrt(t^2 / (N - rank(M) + t^2));
 r_partial{tv,1} = r;
end

cell2mat(r_partial)
tbl.r_partial = r_partial; % concatenate with original table
writetable(tbl, '/home/uni10/nmri/projects/cstier/aging_pipeline/results/global_MEG_age_age2_sex_tiv/r_partial_global_MEG_age_age2_sex_tiv.csv') % save
