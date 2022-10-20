% This scripts loads the statistical results for the age effects on cortical
% thickness and plots the t-values (TFCE-corrected) for significant vertices
% 
% Written by Christina Stier, 2022

% load folders and assign file names
analysis_dir = ('/home/uni10/nmri/projects/cstier/aging_pipeline');
results_dir = '/home/uni10/nmri/projects/cstier/aging_pipeline/results/thickness_age_age2_sex_tiv/';

modality = {'_tfce_tstat_fwep_', '_tfce_tstat_'};
contrast = {'c2', 'c4'}; % choose contrast: here decrease with age (c2) and inverted U-shape with age (c4)

% load settings for plotting  
load([analysis_dir '/conf/suma-all-fsaverage-10.mat'],'suma_all')

type = {'Sign. t-values (tfce, fwep)'};
out = 'cortical_thickness_decrease_';
out2 = {'sign_t_tfce'};

opt=[];
opt.per_hemi=1;
opt.per_cortex=1;
opt.rot=[90 0 ; -90 0];
opt.colormap='summer'; % change to any other color depending on the contrast (hot, viridis, inferno)
opt.colorbar='summer';
opt.scale=1;

% load p and t maps  
filename = fullfile(results_dir, [modality{1} contrast{1} '.mgz']); % choose contrast
p_metric = load_mgh(filename);

t_file = fullfile(results_dir, [modality{2} contrast{1} '.mgz']); % choose contrast
t_metric = load_mgh(t_file);

% select sign. vertices
notsign_val = find(p_metric < 1.3);
t_metric(notsign_val) = 0;

% plot
opt.title= type{1} ;
opt.output= [results_dir, out, out2{1}, '.png'];

hFig = nmri_plot_surface_suma(suma_all, t_metric, opt);
             
 





