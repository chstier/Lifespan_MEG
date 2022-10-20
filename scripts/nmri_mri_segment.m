function [ mri_all ] = nmri_mri_segment( cfg, params )
% [ mri_all ] = nmri_mri_segment( cfg, params )
% will segment (SPM12, new segement) MRI scans for
% use in Fieldtrip (e.g. for headmodeling). Data should be passed in a
% cfg-struct
%
% cfg.t1       =  path to T1 image or Fieltrip mri struct (needed)
% .flair       =  path to FLAIR image or Fieltrip mri struct (optinal)
% .t2          =  path to T2 image or Fieltrip mri struct (optional)
% .write_class = 0: do not write class files, 1: write class files
% .write_norm  = 0: do not write w files, 1: write wfiles
% .dim         = dimensions to use (default [256 256 256])
% 
% It is preferable to have as little image bias as possible, e.g.
% T1.mgz/.nii after Freesurfer / SUMA processing. Or otherwise very "flat"
% / bias corrected scans. SPM12 segement will also be used to reduce bias. 
%
% If multiple channels are used (FLAIR / T2) in addition, to T1, these need
% to be coregisterd FIRST. Otherwise SPM segment will fail.
%

%% Setup of vars and data
if ~exist('cfg','var')
 error('Need a cfg-struct')
end

% check the params struct
if (~exist('params','var') || isempty(params))
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in the current path, or have it in the call ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the parameter file (analysis_params.m)')  
  end
 end
end


if ~isfield(cfg,'t1')
 error('T1 image is always required')
end

if ~isfield(cfg,'write_class')
 cfg.write_class=0;
end
if ~isfield(cfg,'write_seg')
 cfg.write_seg=1;
end
if ~isfield(cfg,'write_norm')
 cfg.write_norm=0;
end
if ~isfield(cfg,'dim')
 cfg.dim=[256 256 256];
end

if ischar(cfg.t1)
 if exist(cfg.t1,'file')
  mri_t1=ft_read_mri(cfg.t1);
 else
  error(['T1 file not found, ' cfg.t1])
 end
elseif isfield(cfg.t1.hdr) && isfield(cfg.t1.anatomy)
 % seems to be a FT MRI struct
 mri_t1=cfg.t1;
else
 error('T1 image is not a char and not an MRI struct')
end


if isfield (cfg,'flair')
 if ischar(cfg.flair)
  if exist(cfg.flair,'file')
   mri_flair=ft_read_mri(cfg.flair);
  else
   error(['FLAIR file not found, ' cfg.flair])
  end
 elseif isfield(cfg.flair.hdr) && isfield(cfg.flair.anatomy)
  % seems to be a FT MRI struct
  mri_flair=cfg.flair;
 else
  error('FLAIR was set, but is not a char and not an MRI struct') 
 end
end

if isfield (cfg,'t2')
 if ischar(cfg.t2)
  if exist(cfg.t2,'file')
   mri_t2=ft_read_mri(cfg.t2);
  else
   error(['T2 file not found, ' cfg.t2])
  end
 elseif isfield(cfg.t2.hdr) && isfield(cfg.t2.anatomy)
  % seems to be a FT MRI struct
  mri_t2=cfg.t2;
 else
  error('T2 was set, but is not a char and not an MRI struct') 
 end
end

% require SPM12 to be present
defaults=nf_request_spm;

% make the paths sets
[T1_pa, T1_fi,ext]=fileparts(mri_t1.hdr.fspec);
T1_fi=[T1_fi ext];
if isempty(T1_pa)
 T1_pa=mri_t1.hdr.pwd;
end

%% Initiate the SPM 12 segementation, if not present
required_classes={'c1','c2','c3','c4','c5','c6','rc1','rc2','y_','iy_','m'};
doseg=0;
% check if we have everything we need
for i=1:length(required_classes)
 if ~exist(fullfile(T1_pa,[required_classes{i} T1_fi]))
  fprintf('Missing class=%s, initiating segmentation\n',required_classes{i})
  doseg=1;
 end
end
if isfield (cfg,'flair')
 [FLAIR_pa, FLAIR_fi, ext]=fileparts(mri_flair.hdr.fspec);
 FLAIR_fi=[FLAIR_fi ext];
 if isempty(FLAIR_pa)
  FLAIR_pa=mri_flair.hdr.pwd;
 end
 if ~exist(fullfile(FLAIR_pa,['m' FLAIR_fi])) && ~exist(fullfile(FLAIR_pa,['mr' FLAIR_fi]))
  fprintf('Missing BIAS corr FLAIR, initiating segmentation\n')
  doseg=1;
 end
end
if isfield (cfg,'t2')
 [T2_pa, T2_fi, ext]=fileparts(mri_t2.hdr.fspec);
 T2_fi=[T2_fi ext];
 if isempty(T2_pa)
  T2_pa=mri_t2.hdr.pwd;
 end
 if ~exist(fullfile(T2_pa,['m' T2_fi])) && ~exist(fullfile(T2_pa,['mr' T2_fi]))
  fprintf('Missing BIAS corr T2, initiating segmentation\n')
  doseg=1;
 end
end

if (doseg==1)
 % run segmentation
 if isfield (cfg,'t2') || isfield (cfg,'flair')
  % multispectral, wants cell array
  files={fullfile(T1_pa,T1_fi)};
  if isfield (cfg,'t2')
   % realign and reslice, if needed
   if ~isequal(mri_t1.transform,mri_t2.transform)
    % need to realign
    nf_coreg_single(4,fullfile(T2_pa,T2_fi),fullfile(T1_pa,T1_fi),T2_pa)
    if ~exist(fullfile(T2_pa,['r' T2_fi]),'file')
     error('Problems while realiging T2 to T1')
    else
     T2_fi=['r' T2_fi];
    end
   end
   files{end+1}=fullfile(T2_pa,T2_fi);
  end
  if isfield (cfg,'flair')
   % realign and reslice, if needed
   if ~isequal(mri_t1.transform,mri_flair.transform)
    % need to realign
    nf_coreg_single(4,fullfile(FLAIR_pa,FLAIR_fi),fullfile(T1_pa,T1_fi),FLAIR_pa)
    if ~exist(fullfile(FLAIR_pa,['r' FLAIR_fi]),'file')
     error('Problems while realiging FLAIR to T1')
    else
     FLAIR_fi=['r' FLAIR_fi];
    end
   end
   files{end+1}=fullfile(FLAIR_pa,FLAIR_fi);
  end
 else
  % T1-only, char array
  files=fullfile(T1_pa,T1_fi);
 end
 % now run it
 nf_batch_segment_new(files,T1_pa,3,1,[ 1 1 1 1 1 1],[ 1 1 0 0 0 0],[ 0 0 0 0 0 0],[ 0 0 0 0 0 0])
else
 fprintf('Segmentation already done...\n')
end


%% Check if DARTEL is also needed
required_classes={'u_rc1'};
doseg=0;
% check if we have everything we need
for i=1:length(required_classes)
 if ~exist(fullfile(T1_pa,[required_classes{i} T1_fi]))
  fprintf('Missing class=%s, initiating DARTEL\n',required_classes{i})
  doseg=1;
 end
end

if (doseg==1)
 % run DARTEL
 clear matlabbatch;
 origdir=fullfile(getenv('NMRI_TOOLS'),'common','DARTEL_templates_CAT12');
 matlabbatch{1}.spm.tools.dartel.warp1.images{1,1} = {fullfile(T1_pa,['rc1' T1_fi])};
 matlabbatch{1}.spm.tools.dartel.warp1.images{1,2} = {fullfile(T1_pa,['rc2' T1_fi])};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(origdir,'Template_1.nii')};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(origdir,'Template_2.nii')};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(origdir,'Template_3.nii')};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(origdir,'Template_4.nii')};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(origdir,'Template_5.nii')};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(origdir,'Template_6.nii')};
 matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
 matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;
 spm_jobman('run',matlabbatch)
else
 fprintf('DARTEL already done...\n')
end

%% Start with Scalp / other shell
% read BIAS-corr T1
mri_all=ft_read_mri(fullfile(T1_pa,['m' T1_fi]));
if isfield(cfg,'dim') && any(mri_all.dim~=cfg.dim)
 rscfg = [];
 rscfg.dim = cfg.dim;
 mri_all = ft_volumereslice(rscfg, mri_all);
end

fprintf('Generating a scalp mask from Freesurfer processed T1...\n')

% get percentile for non-0
msk=mri_all.anatomy;
%msk(msk==0)=[];
% we have a Freesurfer proc T1, so assume 30 as non-tissue and take 5th percentile 
cutoff=prctile(msk(msk>30),5);
fprintf('T1 cutoff=%0.3f\n',cutoff)
clear msk
% smooth mildly
scalp=nf_volumesmooth(mri_all.anatomy,2);
scalp=scalp>cutoff;

% add the c5 class / tissue as most outside class
mri_class=ft_read_mri(fullfile(T1_pa,['c5' T1_fi]));
if isfield(cfg,'dim') && any(mri_class.dim~=cfg.dim)
 rscfg = [];
 rscfg.dim = cfg.dim;
 mri_class = ft_volumereslice(rscfg, mri_class);
end
mri_class.anatomy=nf_volumesmooth(mri_class.anatomy,1);
if ~isfield(params,'c5cutoff') || ~isnumeric(params.c5cutoff)
 c5cutoff=0.8;
else
 c5cutoff=params.c5cutoff;
end


fprintf('Adding the c5 class with prob. cutoff=%0.3f\n',c5cutoff)

scalp=scalp|(mri_class.anatomy>c5cutoff); % but be somewhat strict
% 02-2020: DV suggested a higher threshold here (0.99)
% however this causes issues if there is no information  to due missing
% MRI/field of view. Going back to 0.8 from 09-2020, but allowing to
% configure this in params
% 


% and the the c4 class / skull for good measure
fprintf('Adding the c4 class with prob. cutoff=%0.3f\n',0.8)
mri_class=ft_read_mri(fullfile(T1_pa,['c4' T1_fi]));
if isfield(cfg,'dim') && any(mri_class.dim~=cfg.dim)
 rscfg = [];
 rscfg.dim = cfg.dim;
 mri_class = ft_volumereslice(rscfg, mri_class);
end
mri_class.anatomy=nf_volumesmooth(mri_class.anatomy,1);
scalp=scalp|(mri_class.anatomy>0.8); % but be strict 


% fill holes - first pass

fprintf('Filling holes...\n')
a1 = nf_volumefillholes_all(scalp, 1);
a2 = nf_volumefillholes_all(scalp, 2);
a3 = nf_volumefillholes_all(scalp, 3);
scalp = a1 | a2 | a3;
scalp=nf_volumesmooth(scalp,1);
% re-binarize and also keep only the largest cluster
scalp=nf_volumethreshold(scalp,0.5,'Scalp');
% fill again - but somewhat more mild
a1 = nf_volumefillholes(scalp, 1);
a2 = nf_volumefillholes(scalp, 2);
a3 = nf_volumefillholes(scalp, 3);
scalp = a1 | a2 | a3;
clear a1 a2 a3
fprintf('...done with scalp class\n')
mri_all.scalp=scalp;


% Now we work on classes
required_classes={'c1','c2','c3','c4','c5','c6'};
labels={'grey','white','csf','skull','tissue','air'};
smooth={2,2,2,1,2,2};
weight={1,1,1,1.5,1,1};


% check if we have everything we need
max_v=zeros(mri_all.dim);
for i=1:length(required_classes)
 fprintf('Reading class=%s, prefix=%s\n',labels{i},required_classes{i})
 mri_class=ft_read_mri(fullfile(T1_pa,[required_classes{i} T1_fi]));
 if isfield(cfg,'dim') && any(mri_class.dim~=cfg.dim)
  rscfg = [];
  rscfg.dim = cfg.dim;
  mri_class = ft_volumereslice(rscfg, mri_class);
 end
 msk=mri_class.anatomy;
 if (~isequal(size(msk),mri_all.dim))
  error(['Dimension mismatch for class = ' required_classes{i}]);
 end

 % weigth some classes more
 if weight{i} ~= 1
  msk=msk*weight{i};
 end
 % smooth
 msk=nf_volumesmooth(msk,smooth{i});
 % remove low prob
 msk(msk<0.1)=0;
 % check max
 this_max=msk>max_v;
 max_v(this_max)=msk(this_max);
 % now remove from old classes
 for ii=1:(i-1)
  mri_all.(labels{ii})(this_max)=0;
 end 
 % contain in scalp
 this_max(scalp==0)=0;
 mri_all.(labels{i})=this_max>0;
end

% make a brain mask (in FT logic, including CSF...)
fprintf('Creating brain mask from c1/c2/c3...\n')
brain=mri_all.grey | mri_all.white | mri_all.csf; 
a1 = nf_volumefillholes(brain, 1);
a2 = nf_volumefillholes(brain, 2);
a3 = nf_volumefillholes(brain, 3);
fbrain = a1 | a2 | a3;
fbrain=nf_volumesmooth(fbrain,1);
fbrain=fbrain>0.3;
a1 = nf_volumefillholes(fbrain, 1);
a2 = nf_volumefillholes(fbrain, 2);
a3 = nf_volumefillholes(fbrain, 3);
fbrain = a1 | a2 | a3;
clear a1 a2 a3
% set the difference to be tissue (not skull)
mri_all.brain=fbrain;
diff_brain=(fbrain-brain)>0;
mri_all.tissue(diff_brain)=1;
mri_all.skull(diff_brain)=0;

% now make sure that the skull is closed (foramen magnum needs to be re-opened later)
braindil = imdilate(fbrain>0, strel_bol(2));
skullmask = braindil & ~fbrain;
% set this to skull
mri_all.skull(skullmask)=1;
for i=1:length(labels)
 % remove from all  other classes
 if ~strcmp(labels{i},'skull')
  mri_all.(labels{i})(skullmask)=0;
 end
end

clear brain fbrain diff_brain skullmask braindil

% now we apply deformation to T1 (with generic bounding box)
clear matlabbatch
% NOTE: this is based on the iy_ of the UNIFIED Segmentation! Not DARTEL
matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(T1_pa,['iy_' T1_fi])};
matlabbatch{1}.spm.util.defs.out{1}.push.fnames={};
if (cfg.write_norm==1)
 for i=1:length(required_classes)
  matlabbatch{1}.spm.util.defs.out{1}.push.fnames {end+1,1} = fullfile(T1_pa,[required_classes{i} T1_fi]);
 end
end
matlabbatch{1}.spm.util.defs.out{1}.push.fnames {end+1,1}= fullfile(T1_pa,['m' T1_fi]);
matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
matlabbatch{1}.spm.util.defs.out{1}.push.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb = [-120 -140 -150
                                                         120 120 120];% bounding box to be used by SPM deformations - empirically defined
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox = [1 1 1];
matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.push.prefix = 'w';

spm_jobman('run',matlabbatch);


% now load the MNI norm T1 and make an MNI mask
fprintf('Preparing MNI mask\n')
mni_t1_vol=spm_vol(fullfile(T1_pa,['wm' T1_fi]));
mni_t1_vol.fname=fullfile(T1_pa,'MNImask.nii');
spm_write_vol(mni_t1_vol,ones(mni_t1_vol.dim));
clear mni_t1_vol

% now deform back to original
clear matlabbatch
matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(T1_pa,['iy_' T1_fi])};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(T1_pa,'MNImask.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'i_w';

spm_jobman('run',matlabbatch);

% now load the mask
mni_msk_vol=spm_vol(fullfile(T1_pa,'i_wMNImask.nii'));
mni_msk=spm_read_vols(mni_msk_vol);
mni_msk=mni_msk>0.5;

% now the same jiggle for a tighter skull mask
% now we apply deformation to T1 (with generic bounding box)
clear matlabbatch
matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(T1_pa,['iy_' T1_fi])};
matlabbatch{1}.spm.util.defs.out{1}.push.fnames {1}= fullfile(T1_pa,['m' T1_fi]);
matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
matlabbatch{1}.spm.util.defs.out{1}.push.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb = [NaN NaN NaN; NaN NaN NaN];% bounding box of the priors
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox = [1 1 1];
matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.push.prefix = 'priors_w';

spm_jobman('run',matlabbatch);


% now load the MNI norm T1 and make an MNI mask
fprintf('Preparing MNI mask (tight)\n')
mni_t1_vol=spm_vol(fullfile(T1_pa,['priors_wm' T1_fi]));
mni_t1_vol.fname=fullfile(T1_pa,'priorsMNImask.nii');
spm_write_vol(mni_t1_vol,ones(mni_t1_vol.dim));
clear mni_t1_vol

% now deform back to original
clear matlabbatch
matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(T1_pa,['iy_' T1_fi])};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(T1_pa,'priorsMNImask.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'i_w';

spm_jobman('run',matlabbatch);

% now load the mask
mni_msk_vol=spm_vol(fullfile(T1_pa,'i_wpriorsMNImask.nii'));
mni_msk_priors=spm_read_vols(mni_msk_vol);
mni_msk_priors=mni_msk_priors>0.5;


% now apply to skull and make overlap as tissue
diff_skull=(mri_all.skull&(~mni_msk_priors));
mri_all.tissue(diff_skull)=1;
mri_all.skull(diff_skull)=0;

% now make final masking (MNI derived) and save
labels={'grey','white','csf','skull','tissue','brain','scalp','air'};
% use SPM to safe... Fieldtrip tries to be "smart" with the fields...
t1_vol=spm_vol(fullfile(T1_pa,T1_fi));
for i=1:length(labels)
 fprintf('Finalizing class = %s\n',labels{i})
 t1_vol.fname=fullfile(T1_pa,[labels{i} '.nii']);
 if ~strcmp(labels{i},'air')
  mri_all.(labels{i})(~mni_msk)=0;
 else
  % make everything non-scalp as air
  mri_all.air(mri_all.scalp==0)=1;
 end
 if (cfg.write_class==1)
  data      = mri_all.(labels{i})>0; 
  fprintf('Writing class = %s\n',labels{i})
  spm_write_vol(t1_vol,data);
 end
end

% write out indexed segmentation and pseudo-T1, if requested
if (cfg.write_seg==1)
 labels={'grey','white','csf','skull','tissue'};
 class_idx=[1 2 3 4 5];
 class_val=[50 100 30 10 200 ];
 mri_idx=zeros(mri_all.dim);
 mri_T1=zeros(mri_all.dim);
 for i=1:length(labels)
  mri_idx(mri_all.(labels{i}))=class_idx(i);
  mri_T1(mri_all.(labels{i}))=class_val(i);
 end
 fprintf('Writing indexed segmentation\n')
 t1_vol.fname=fullfile(T1_pa,'segmentation_idx.nii');
 spm_write_vol(t1_vol,mri_idx);
 fprintf('Writing pseudo-T1\n')
 t1_vol.fname=fullfile(T1_pa,'segmentation_pseudoT1.nii');
 spm_write_vol(t1_vol,mri_T1);
end


% cleanup
system(['rm -f ' fullfile(T1_pa,'MNImask.nii') ' ' fullfile(T1_pa,'i_wMNImask.nii') ' ' fullfile(T1_pa,'i_wpriorsMNImask.nii') ' ' fullfile(T1_pa,'priorsMNImask.nii') ' ' fullfile(T1_pa,['priors_wm' T1_fi])]);


end

