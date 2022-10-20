function [cu_suma_all] = cs_get_thickness_suma_all_sm(subject, opt, suma_folder, mod)
%% This function loads cortical thickness values for each subject and concatenates
% hemispheres. In addition, it plots the map and stores in subject-folder
% Written and modified by Christina Stier, 2020/22

ld = opt.ld; % choose density factor
suma_dir = fullfile(suma_folder, [subject.id mod], 'SUMA');

% load thickness
cfile_LH = fullfile(suma_dir, opt.left_file); 
cfile_RH = fullfile(suma_dir, opt.right_file);

dest_LH = afni_niml_read(cfile_LH); % based on AFNI matlab scripts, should be in utilities
dest_RH = afni_niml_read(cfile_RH); % based on AFNI matlab scripts, should be in utilities

% get thickness-values for 1002 vertices for each hemisphere
th_left = dest_LH{1,1}.nodes{1,1}.data;
th_right = dest_RH{1,1}.nodes{1,1}.data;

% concat thickness of both hemispheres
th_cortex = vertcat(th_left, th_right);

% add zeros for subcortical suma points (can be present in functional data)
no_data = zeros(str2num(opt.subc), 1); % ld 10 = 334; ld 40 = 5390
cu_suma_all = vertcat(th_cortex, no_data);

%% plot and save in QC-folder of the subject

if ~exist(fullfile(subject.QCdir, ['thickness_' ld '_' subject.id '_sigm_' opt.sigma '_fwhm_' opt.fwhm '.png']), 'file')
 disp('plotting thickness...')
 load([pwd '/conf/suma-all-fsaverage-' ld '.mat'],'suma_all')
 target_dir = subject.QCdir;

 opt2=[];
 opt2.per_hemi=1;
 opt2.per_cortex=1;
 opt2.rot=[90 0 ; -90 0];
 opt2.thresh=0;
 opt2.clim=[0 4];
 opt2.colormap='hot';
 opt2.colorbar='hot';
 opt2.scale=1;
 opt2.title=fullfile(['Cortical thickness ' subject.id ' smoothed' ]);
 opt2.output=fullfile(target_dir, ['thickness_' ld '_' subject.id '_sigm_' opt.sigma '_fwhm_' opt.fwhm '.png']);

 hFig = nmri_plot_surface_suma_test(suma_all, cu_suma_all, opt2);
else
 disp('thickness-plot already present')
end




 
end 