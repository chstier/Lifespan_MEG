function [cu_suma_all] = cs_get_thickness_suma_all_nosm(subject, opt)

suma_dir = fullfile('/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/', [subject.id '_anat_T2_T1'], 'SUMA');
ld = opt.ld; % choose density factor

% load thickness
cfile_LH = fullfile(suma_dir, opt.left_file); 
cfile_RH = fullfile(suma_dir, opt.right_file);

% % try with _cap if not found
% if ~exist(cfile_LH,'file')
%  warning(['SUMA thickness not found = ' cfile_LH])
%  suma_dir = fullfile('/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/', [subject.id '_align_big_bound_cap'], 'SUMA');
%  cfile_LH = fullfile(suma_dir, opt.left_file);
%  cfile_RH = fullfile(suma_dir, opt.right_file);
% end
% 
% % try with nu_corr_cap if not found
% if ~exist(cfile_LH,'file')
%  warning(['SUMA thickness not found = ' cfile_LH])
%  suma_dir = fullfile('/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/', [subject.id '_anat_T2_T1_cap'], 'SUMA');
%  cfile_LH = fullfile(suma_dir, opt.left_file);
%  cfile_RH = fullfile(suma_dir, opt.right_file);
% end
% 
% % delete letter for timepoint and try again; sometimes there is a mismatch
% % between functional and structural labels
% if ~exist(cfile_RH,'file')
%  warning(['SUMA thickness not found = ' cfile_RH])
%  subj_name = subject.id(1:end-1)
%  fs_folders = dir(fullfile(['/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/' subj_name '*'])); % look for similar ID 
%  t1_folder = fs_folders(1).name(1:9); % take new ID
%  suma_dir = fullfile('/home/uni10/nmri/projects/cstier/aging_freesurfer/6.0.0/', [t1_folder '_anat_T2_T1'], 'SUMA');
%  cfile_LH = fullfile(suma_dir, opt.left_file);
%  cfile_RH = fullfile(suma_dir, opt.right_file);
% end
% 
% % still not found?
% if ~exist(cfile_RH,'file')
%  error(['SUMA thickness still not found = ' cfile_RH])
% end

dest_LH = afni_niml_read(cfile_LH); % based on AFNI matlab scripts, should be in utilities
dest_RH = afni_niml_read(cfile_RH); % based on AFNI matlab scripts, should be in utilities

% get thickness-values for 1002 vertices for each hemisphere
th_left = dest_LH{1,1}.nodes{1,1}.data;
th_right = dest_RH{1,1}.nodes{1,1}.data;


% concat thickness of both hemispheres
th_cortex = vertcat(th_left, th_right);

% add zeros for subcortical suma points (present in functional data)
no_data = zeros(str2num(opt.subc), 1); % ld 10 = 334; ld 40 = 5390
cu_suma_all = vertcat(th_cortex, no_data);

%% plot and save in QC-folder of the subject

if ~exist(fullfile(subject.QCdir, ['thickness_' ld '_nosm_' subject.id '.png']), 'file')
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
 opt2.title=fullfile(['Cortical thickness ' subject.id]);
 opt2.output=fullfile(target_dir, ['thickness_' ld '_nosm_' subject.id '.png']);

 hFig = nmri_plot_surface_suma_test(suma_all, cu_suma_all, opt2);
else
 disp('thickness-plot already present')
end




 
end 