function [ subject ] = nmri_align_mri(subject, params)
%[ subject ] = nmri_align_mri(subject, params)
%  
% This function will do (interactive) MRI alignement to CTF 
% and also define the SUMA dir
%
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file

% written by NF 11/2016


% Freesurfer storage dir
if ~isempty(getenv('SUBJECTS_DIR'))
 fsdir = getenv('SUBJECTS_DIR');
else
 % default to NMRI/PATLAN behaviour
 fsdir = '/data/freesurfer/6.0.0';
end


% check the LR-marker image
if ~isempty(getenv('NMRI_TOOLS'))
 LR_marker_image = fullfile(getenv('NMRI_TOOLS'),'common','freesurfer-SUMA-L-R-mask');
else
 LR_marker_image='/tools/common/freesurfer-SUMA-L-R-mask';
end

% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps


%% check the SUMA dir - if we do have all MRIs/SUMA yet
if (~isfield(subject,'suma_dir') && ( ~exist(subject.mri_ctf,'file') || ~exist(subject.mri_seg,'file') || ~exist(subject.suma_surface,'file')))
 % search for matching fressurfer IDs
 res=dir(fullfile(fsdir,[ subject.id '*']));
 possible={};
 for i=1:length(res)
  % check SUMA dir
  if (exist(fullfile(fsdir,res(i).name,'SUMA','std.40.lh.white.gii'),'file'))
   %seems processed
   possible{end+1}=res(i).name;
  else
   disp(fullfile(fsdir,res(i).name,'SUMA','std.40.lh.white.gii'))
  end
 end
 if (length(possible)==1)
  % only one, ask user if happy
  qcfg=[];
  qcfg.question={['Found one match=' possible{1}],'Are you happy with this choice?','Hit ''No'' to manually select.'};
  qcfg.title='SUMA found';
  qcfg.options={'Yes','No'};
  qcfg.default={'Yes'};
  qcfg.mode='buttons';
  button=nf_uidialog(qcfg);
 
  if strcmp(button,'Yes')
   % then take it
   subject.fsid=possible{1};
   subject.suma_dir=fullfile(fsdir,subject.fsid,'SUMA');
   subject.fsdir=fullfile(fsdir,subject.fsid);
  end
 end
 if (length(possible)>1)
  % more than one, ask user to pick
  qcfg=[];
  qcfg.question={['Found >1 possibility, please pick one (Subject=' subject.id  ')']};
  qcfg.title='SUMA selection';
  qcfg.options=possible;
  qcfg.mode='popup';
  manselect=nf_uidialog(qcfg);
  %manselect=listdlg('Name',,'SelectionMode','single','ListSize',[300 (50+(length(possible)*10))],'ListString',possible);
  if (~isempty(manselect))
   subject.fsid=manselect;
   subject.suma_dir=fullfile(fsdir,subject.fsid,'SUMA');
   subject.fsdir=fullfile(fsdir,subject.fsid);
  end
 end
 % if nothing, ask user
 while (~isfield(subject,'suma_dir'))
  manselect=uigetdir(fsdir,['Nothing found - Select subject dir manually (Subject=' subject.id  ')']);
  if (manselect==0)
   error('Nothing selecting, exiting')
  end
  if (exist(fullfile(manselect,'SUMA','std.40.lh.white.gii'),'file'))
   subject.suma_dir=fullfile(manselect,'SUMA');
   subject.fsdir=manselect;
   [pa,fi,ext]=fileparts(manselect);
   subject.fsid=[fi ext];
  else
   disp(['Could not find processed SUMA files for this subject = ' fullfile(fsdir,manselect,'SUMA','std.40.lh.white.gii') ', pick another one'])
  end
 end 
 disp(sprintf('Selected SUMA dir = %s',subject.suma_dir))
end


% make our output paths and dir

if (~isfield(subject,'mri_dir'))
 subject.mri_dir=fullfile(subject.analysis_dir,subject.id,'mri');
end
if (~exist(subject.mri_dir,'dir'))
 mkdir(subject.mri_dir)
end

if (~isfield(subject,'mri_ctf'))
 subject.mri_ctf=fullfile(subject.mri_dir,['mri_ctf_' subject.id '.mat']);
end
if (~isfield(subject,'mri_T1'))
 subject.mri_T1=fullfile(subject.mri_dir,['mri_T1_' subject.id '.nii']);
end


%% Read MRI

% Load T1 and as to manually mark LPA - RPA - nasion and Z-dir
if (exist(subject.mri_ctf,'file'))
 disp('Found CTF-aligend MRI - nothing to do...')
else
 % not present - Process the Freesurfer / SUMA T1
 disp('CTF-aligend MRI not found - processing')
 
 %t1_file=fullfile(subject.suma_dir,'T1.nii');
 t1=fullfile(subject.fsdir,'mri','nu.mgz'); % we now take the nu _cor T1
 t1_file=subject.mri_T1;
 t1L_file=fullfile(subject.mri_dir,['mri_T1LR' subject.id '.nii.gz']);
 if ~exist(t1_file,'file') || ~exist(t1L_file,'file')
  if ~exist(t1,'file')
   error(['Could not find Freesurfer image T1=' t1])
  else
   [status cmdout]=system(['mri_convert ' t1 ' ' subject.mri_T1]);
   if (status~=0 || ~exist(subject.mri_T1,'file'))
    cmdout
    error('Could not copy Freesurfer T1 image, likely permission or filesystem problem')
   end
  
   % add the LR-marker via FSL - this should be orientation safe for NIFTI
   % conforming datasets
   [status cmdout]=system(['fslmaths ' subject.mri_T1 ' -add ' LR_marker_image ' ' t1L_file]);
   if (status==0 && exist(t1L_file,'file'))
    t1_file=t1L_file;   
   else
    warning('Could not add the LR-side marker, something may be wrong with your T1 file or the FSL installation')
    cmdout
   end
  end
 end
 % check outside copy loop
 if exist(t1L_file,'file')
  t1_file=t1L_file;   
 end
 mri = ft_read_mri(t1_file,'dataformat','nifti');
  
 cfg = [];
 if isfield(params,'MRI') && isfield(params.MRI,'dim')
  cfg.dim = params.MRI.dim;
 else
  cfg.dim  = [256 256 256];
 end
 mri_rs = ft_volumereslice(cfg, mri); % applies the read in transformation matrix
 mri_rs = ft_convert_units(mri_rs, 'cm'); % conform to MEG units
  
 % realign to the CTF sensor space using the fiducials
 cfg = [];
 cfg.method  = 'interactive';
 mri_ctf = ft_volumerealign(cfg, mri_rs);

  
 % check if we have resampled and have fiducials
 if (isfield(mri_ctf,'transformorig') && isfield(mri_ctf,'cfg') && isfield(mri_ctf.cfg,'fiducial') && ~any(isnan(mri_ctf.cfg.fiducial.lpa)) && ~any(isnan(mri_ctf.cfg.fiducial.rpa)) && ~any(isnan(mri_ctf.cfg.fiducial.nas)))
  save(subject.mri_ctf, 'mri_ctf','mri');
  % save the suma_dir in dws_filt_dataset
  if isfield(subject,'dws_filt_dataset') && exist(subject.dws_filt_dataset,'file')
   subject_old=load(subject.dws_filt_dataset,'subject');
   if ~isequal(subject_old.subject,subject)
    update_dirs(subject_old.subject,subject,subject.dws_filt_dataset)
   end 
  end
  % save the suma_dir in clean_dataset
  if isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
   subject_old=load(subject.clean_dataset,'subject');
   if ~isequal(subject_old.subject,subject)
    update_dirs(subject_old.subject,subject,subject.clean_dataset)
   end 
  end
  % save the suma_dir in cleanICA_dataset
  if isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')
   subject_old=load(subject.cleanICA_dataset,'subject');
   if ~isequal(subject_old.subject,subject)
    update_dirs(subject_old.subject,subject,subject.cleanICA_dataset)
   end 
  end
  
  % if we got this far, place a stamp for completed MRI-alignement
  subject=nmri_stamp_subject(subject,'MRI_ctf',params);
 else
  error('It seems you have not specified the markers (LPA,RPA,nasion are minimum), no transformation was done')
 end
end


end

function update_dirs(subject,subject_updated,fname)
 % update the relevant fields, leave rest intact
 items=fieldnames(subject_updated);
 for i=1:length(items)
  if length(items{i})>=2 && (strcmp(items{i}(1:2),'fs'))
   subject.(items{i})=subject_updated.(items{i});
  end
  if length(items{i})>=4 && (strcmp(items{i}(1:4),'suma'))
   subject.(items{i})=subject_updated.(items{i});
  end
 end
 fprintf('Updating SUMA and Freesurfer info in %s\n',fname)
 save(fname,'subject','-append');
end

