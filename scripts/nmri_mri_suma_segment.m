function [ subject ] = nmri_mri_suma_segment(subject, params)
% [ subject ] = nmri_mri_suma_segment(subject, params)
%  
% This function will deal with reading the MRI info 
% per default will read from SUMA, other options may be integrated later 
% will assume that the subject.id = basename
% otherwise you should specify subject.suma_dir
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file

% written by NF 11/2016


% Freedurder storage dir
fsdir = getenv('SUBJECTS_DIR');


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% make our output paths and dir

if (~isfield(subject,'mri_dir'))
 subject.mri_dir=fullfile(subject.analysis_dir,subject.id,'mri');
end
if (~exist(subject.mri_dir,'dir'))
 mkdir(subject.mri_dir)
end

if (~isfield(subject,'proc_dir'))
 subject.proc_dir=fullfile(subject.analysis_dir,subject.id,'processed');
end
if (~exist(subject.proc_dir,'dir'))
 mkdir(subject.proc_dir)
end

% These are subject specific, so could be re-used for all exams
if (~isfield(subject,'mri_ctf'))
 subject.mri_ctf=fullfile(subject.mri_dir,['mri_ctf_' subject.id '.mat']);
end
if (~isfield(subject,'mri_seg'))
 subject.mri_seg=fullfile(subject.mri_dir,['mri_seg_' subject.id '.mat']);
end
if (~isfield(subject,'mri_T1'))
 subject.mri_T1=fullfile(subject.mri_dir,['mri_T1_' subject.id '.nii']);
end
if (~isfield(subject,'suma_surface'))
 subject.suma_surface=fullfile(subject.mri_dir,['suma_surface_' subject.id '.mat']);
end

% This is subject and modality specific
if (~isfield(subject,'hdm_lead'))
 subject.hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.dtype '_' params.headmodel '.mat']);
end


%% Read MRI - interactive part

% Load T1 and as to manually mark LPA - RPA - nasion and Z-dir
if (exist(subject.mri_ctf,'file'))
 disp('Found CTF-aligend MRI - loading')
 load(subject.mri_ctf,'mri_ctf','mri')
else
 % not present - Call interactive subscript
 subject=nmri_align_mri(subject);
 if (exist(subject.mri_ctf,'file'))
  load(subject.mri_ctf,'mri_ctf','mri')
 else
  error('Problems with getting the CTF-aligned MRI')
 end
end
  
%% Now comes the processing - non interactive - MRI segmentation is first

% require SPM12 to be present
defaults=nf_request_spm;

% check the MRI segmentation
if (exist(subject.mri_seg,'file'))
 disp('Found segmented MRI - loading')
 load(subject.mri_seg,'mri_seg')
else
 % not present - do the SPM segemenation
 cfg=[];
 cfg.t1=subject.mri_T1;
 cfg.write_norm=1;
 cfg.write_class=0; % rather safe as mat-file later, more compact
 if (isfield(params,'use_FLAIR') && params.use_FLAIR==1)
  subject.mri_FLAIR=fullfile(subject.mri_dir,['mri_FLAIR_' subject.id '.nii']);
  if  ~exist(subject.mri_FLAIR,'file')
   % we take the orignal FLAIR sent to T1-Freesrufer space
   flair=fullfile(subject.fsdir,'mri','orig','FLAIRraw.mgz'); % for FS 6.0.0 
   flair_reg=fullfile(subject.fsdir,'mri','transforms','FLAIRraw.lta'); % for FS 6.0.0   
   if exist(flair,'file') && exist(flair_reg,'file')
    [status cmdout]=system(['mri_vol2vol --mov ' flair ' --reg ' flair_reg ' --fstarg --o ' subject.mri_FLAIR]);
    if (status==0 && exist(subject.mri_FLAIR,'file'))
     cfg.flair=subject.mri_FLAIR;
    else
     cmdout
     error('Could not copy FLAIR image, likely permission or filesystem problem')
    end
   elseif (isfield(params,'use_FLAIR_permissive') && params.use_FLAIR_permissive==1)  
    warning('A FLAIR file was requested, but no FLAIR was found. Permissive mode, will continue without FLAIR')
   else
    error('A FLAIR file was requested, but no FLAIR was found. Change to permissive mode, if you still want to proceed.')
   end
  else
   cfg.flair=subject.mri_FLAIR;
  end
 end
 % T2 
 if (isfield(params,'use_T2') && params.use_T2==1)
  subject.mri_T2=fullfile(subject.mri_dir,['mri_T2_' subject.id '.nii']);
  if ~exist(subject.mri_T2,'file')
   % we take the orignal T2 sent to T1-Freesrufer space
   t2=fullfile(subject.fsdir,'mri','orig','T2raw.mgz'); % for FS 6.0.0 
   t2_reg=fullfile(subject.fsdir,'mri','transforms','T2raw.lta'); % for FS 6.0.0   
   if exist(t2,'file') && exist(t2_reg,'file')
    [status cmdout]=system(['mri_vol2vol --mov ' t2 ' --reg ' t2_reg ' --fstarg --o ' subject.mri_T2]);
    if (status==0 && exist(subject.mri_T2,'file'))
     cfg.t2=subject.mri_T2;
    else
     cmdout
     error('Could not copy T2 image, likely permission or filesystem problem')
    end
   elseif (isfield(params,'use_T2_permissive') && params.use_T2_permissive==1)  
    warning('A T2 file was requested, but no T2 was found. Permissive mode, will continue without T2')
   else
    error('A T2 file was requested, but no T2 was found. Change to permissive mode, if you still want to proceed.')
   end
  else
   cfg.t2=subject.mri_T2;
  end
 end
 if isfield(params,'MRI') && isfield(params.MRI,'dim') 
  cfg.dim=params.MRI.dim;
 end

 % note, this MRI_seg is not fully fieldtrip compatible (air/tissue/scalp).
 % This needs to be sorted out later
 [ mri_seg ] = nmri_mri_segment( cfg );
 save(subject.mri_seg, 'mri_seg');
end

%% now comes the SUMA / ASEG mapping
if (exist(subject.suma_surface,'file'))
 disp('Found SUMA surface - loading')
 load(subject.suma_surface,'suma_all');
else
 % process SUMA

 % load raw SUMA
 suma_all=nf_read_suma(subject.suma_dir,params.SUMA_ld,params.SUMA_surf);
 % load ASEG, if requestes
 if (isfield(params,'ASEG_use') && params.ASEG_use==1)
  if (isfield(params,'ASEG_individual') && params.ASEG_individual==1)
   % this will run an indicidual aseg ROI generation
   aseg_all=nmri_read_aseg(subject.fsdir,params.ASEG_conf);
   % now apply initial MRI transform - then we should be in same space as
   % SUMA
   aseg_all = ft_transform_geometry(vox2ras_0to1(mri.transform), aseg_all);   
  elseif isfield(params,'ASEG_MNI_unified') && params.ASEG_MNI_unified==1
   % based on unified segmenation 
   % load the fsaverge ROIs for this ld
   r_file=fullfile(subject.analysis_dir,'conf',['freesurfer-fsaverage-ROIs-aseg-connectome-' params.SUMA_ld  '.mat']);
   fs_mni=fullfile(subject.analysis_dir,'conf',['freesurfer-fsaverage-ROIs-aseg-connectome-MNI-' params.SUMA_ld  '.gii']);
   if (exist(r_file,'file') && exist(fs_mni,'file')) 
    load(r_file,'surf_avg_ras')
   else
    error('Could not find the pre-computed fsaverage ROIs of this subject (GIFTI and .mat file needed)')
    % if you need to generate for this LD:
    % surf_avg=nf_read_aseg('/data/freesurfer/5.3.0/fsaverage',conf_file)
    % surf_avg_ras=ft_transform_geometry(surf_avg.trans_aseg,surf_avg)
    % save('freesurfer-fsaverage-ROIs-aseg-connectome-YOURLD.mat','surf_avg','surf_avg_ras')
    % ft_write_headshape('freesurfer-fsaverage-ROIs-aseg-connectome-YOURLD.gii',surf_avg_ras,'format','gifti')
    % then apply SPM12 deformation iy_ (to MNI)
   end
   % now apply the individual MNI-native transform
   [pa, fi, ext]=fileparts(subject.mri_T1);
   y_file=fullfile(pa,['y_' fi ext]);
   clear matlabbatch
   if (~exist(y_file,'file'))
    error('Could not find the y_ (MNI_to_native) transformation')
   end
   matlabbatch{1}.spm.util.defs.comp{1}.def = {y_file};
   matlabbatch{1}.spm.util.defs.out{1}.surf.surface = {fs_mni};
   matlabbatch{1}.spm.util.defs.out{1}.surf.savedir.saveusr = {subject.mri_dir};
   spm_jobman('run',matlabbatch);
   fs_to_nat=fullfile(subject.mri_dir,['freesurfer-fsaverage-ROIs-aseg-connectome-MNI-' params.SUMA_ld  '_warped.gii']);
   if (~exist(fs_to_nat,'file'))
    error('Could not generate the warped surface')
   end
 
   aseg_all=surf_avg_ras;
   temp=ft_read_headshape(fs_to_nat);
   aseg_all.pos=temp.pos; 
   clear surf_avg_ras temp  
   
  else 
   % based on DARTEL - new default as of 20170412
   % load the fsaverge ROIs for this ld
   r_file=fullfile(subject.analysis_dir,'conf',['freesurfer-fsaverage-ROIs-aseg-connectome-' params.SUMA_ld  '.mat']);
   fs_mni=fullfile(subject.analysis_dir,'conf',['freesurfer-fsaverage-ROIs-aseg-connectome-DARTEL-' params.SUMA_ld  '.gii']);
   if (exist(r_file,'file') && exist(fs_mni,'file')) 
    load(r_file,'surf_avg_ras')
   else
    error('Could not find the pre-computed fsaverage ROIs of this subject (GIFTI and .mat file needed)')
   end
   % now apply the individual MNI-native transform
   [pa, fi, ext]=fileparts(subject.mri_T1);
   u_file=fullfile(pa,['u_rc1' fi ext]);
   if (~exist(u_file,'file'))
    error('Could not find the u_rc1 (DARTEL) transformation')
   end
   clear matlabbatch
   matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {u_file};
   matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
   matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
   matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};
   matlabbatch{1}.spm.util.defs.out{1}.surf.surface = {fs_mni};
   matlabbatch{1}.spm.util.defs.out{1}.surf.savedir.saveusr = {subject.mri_dir};
   spm_jobman('run',matlabbatch);
   fs_to_nat=fullfile(subject.mri_dir,['freesurfer-fsaverage-ROIs-aseg-connectome-DARTEL-' params.SUMA_ld  '_warped.gii']);
   if (~exist(fs_to_nat,'file'))
    error('Could not generate the warped surface')
   end
 
   aseg_all=surf_avg_ras;
   temp=ft_read_headshape(fs_to_nat);
   aseg_all.pos=temp.pos; 
   clear surf_avg_ras temp  
  end
  
  % add some curv - if we have the toolbox
  if (exist('GetCurvatures','file')==2)
   % found toolbox
   cs.vertices=aseg_all.pos;
   cs.faces=aseg_all.tri;
   curv=GetCurvatures(cs,0);
   aseg_all.curv=(mean(curv,1))';
  else
   % add empty curv 
   aseg_all.curv=zeros(size(aseg_all.pos,1),1);
  end
  if (isfield(aseg_all,'trans_aseg'))
   % remove transformation, to avoid warning
   aseg_all=rmfield(aseg_all,'trans_aseg');
  end
  % now concatenateand 
  suma_all=nf_concat_surf(suma_all,aseg_all);
 end
 
 % also save as GIFTI surface before CTF-transform
 ft_write_headshape(fullfile(subject.mri_dir,['suma_all_' subject.id '.gii']),suma_all,'format','gifti');
 
 % convert to cm
 suma_all = ft_convert_units(suma_all, 'cm');
 
 % note the transformation difference of manual alignement
 T = mri_ctf.transform*inv(mri_ctf.transformorig); 
 % and apply to have in CTF
 suma_all = ft_transform_geometry(T, suma_all);
 suma_all.coordsys='ctf';
 % and save SUMA surface
 save(subject.suma_surface,'suma_all');

end


%% if we got this far, place a stamp for completed SUMA
subject=nmri_stamp_subject(subject,'MRI_seg_SUMA',params);
 

disp('...MRI processing done')
disp('If you want you can now move this subjects to another server (e.g. HIH server) and continue there'); 
 
end

