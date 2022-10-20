function [ suma ] = nf_read_suma( suma_dir, ld, surf_class)
%[ suma ] = nf_read_suma( suma_dir, ld, surf_class)
%   Will read in the SUMA surface of a specific ld from the specified dir
%   will return a strcut with
% .pos = vertices
% .tri = faces
% .curv = curvature
% .annot = annotations
% .annot_key = annotations keys
%
% suma_dir  = path to the SUMA folder (usually <freesurfer-subject>/SUMA)
% ld        = SUMA ld factor
% surf_class= type of surface to use (pial, white, smoothwm, inflated,
%             sphere, mid-cortex)
%             note: mid-cortex is auto-generated to be halfway between pial
%             and white
%
% requires AFNI functions in the path (fiedtrip should take care of that,
% otherwise do yourself)
%
% NF 11/2016

if ~exist('suma_dir','var')
 error('Need SUMA dir')
end

if ~exist('ld','var')
 error('Need SUMA ld factor')
end

if ~exist('surf_class','var')
 error('Need SUMA surface class (pial, white, smoothwm, inflated, sphere, mid-cortex)')
end


if ~strcmp(surf_class,'pial') && ~strcmp(surf_class,'white') && ~strcmp(surf_class,'smoothwm') && ~strcmp(surf_class,'inflated') && ~strcmp(surf_class,'sphere') && ~strcmp(surf_class,'mid-cortex')
 error('No valid SUMA surface class (pial, white, smoothwm, inflated, sphere, mid-cortex)')
end

if isnumeric(ld)
 % make char
 ld=num2str(ld);
end

if ~exist(suma_dir,'dir')
 error('SUMA dir not found')
end


% load SUMA - will always do since non-interactive

if strcmp(surf_class,'mid-cortex')
 % special case to use the mid between gm and wm
 % recurse to get both
 fprintf('Reading pial...\n')
 [ suma ] = nf_read_suma( suma_dir, ld , 'pial');
 fprintf('Reading white...\n')
 [ suma_white ] = nf_read_suma( suma_dir, ld , 'white');
 fprintf('Making mid-cortex pial...\n')
 % and avarge
 
 suma.pos=(suma.pos+suma_white.pos)/2;
 suma.curv=(suma.curv+suma_white.curv)/2;

else
 % all other cases load from Freesufer/SUMA directly
 file_LH = fullfile(suma_dir,['std.' ld '.lh.' surf_class  '.gii']); % change to .gii in FS 6.0.0
 file_RH = fullfile(suma_dir,['std.' ld '.rh.' surf_class  '.gii']);
 if ~exist(file_LH,'file')
  error(['SUMA surface not found = ' file_LH])
 end
 if ~exist(file_RH,'file')
  error(['SUMA surface not found = ' file_RH])
 end

 [suma_LH.pos suma_LH.tri] = nf_load_surf_flexibel(file_LH);
 [suma_RH.pos suma_RH.tri] = nf_load_surf_flexibel(file_RH);

 % load curv

 cfile_LH = fullfile(suma_dir,['std.' ld '.lh.curv.niml.dset']);
 cfile_RH = fullfile(suma_dir,['std.' ld '.rh.curv.niml.dset']);

 if ~exist(cfile_LH,'file')
  error(['SUMA curv not found = ' cfile_LH])
 end
 if ~exist(cfile_RH,'file')
  error(['SUMA curv not found = ' cfile_RH])
 end

 dest_LH = afni_niml_read(cfile_LH); % based on AFNI matlab scripts, should be in utilities
 dest_RH = afni_niml_read(cfile_RH); % based on AFNI matlab scripts, should be in utilities

 % NOTE: the behaviour of Fieldtrip's afni_niml_read is different - does not give cells 
 if ~iscell(dest_LH)
  dest_LH{1}=dest_LH;
 end
 if ~iscell(dest_RH)
  dest_RH{1}=dest_RH;
 end



 if strcmp(dest_LH{1}.nodes{1}.name,'SPARSE_DATA') && isa(dest_LH{1}.nodes{1}.data,'double') % safety check
  curv_LH = dest_LH{1}.nodes{1}.data; % assumes that curv is in first node, usually the case
 else
  error('Unexpected trouble with reading the curvature file - LH')
 end
 if strcmp(dest_RH{1}.nodes{1}.name,'SPARSE_DATA') && isa(dest_RH{1}.nodes{1}.data,'double') % safety check
  curv_RH = dest_RH{1}.nodes{1}.data; % assumes that curv is in first node, usually the case
 else
  error('Unexpected trouble with reading the curvature file - RH')
 end

 % load labels

 afile_LH = fullfile(suma_dir,['std.' ld '.lh.aparc.a2009s.annot.niml.dset']);
 afile_RH = fullfile(suma_dir,['std.' ld '.rh.aparc.a2009s.annot.niml.dset']);

 if ~exist(afile_LH,'file')
  error(['SUMA annotations not found = ' afile_LH])
 end
 if ~exist(afile_RH,'file')
  error(['SUMA annotations not found = ' afile_RH])
 end

 dannot_LH = afni_niml_read(afile_LH); % based on AFNI matlab scripts, should be in utilities
 dannot_RH = afni_niml_read(afile_RH); % based on AFNI matlab scripts, should be in utilities

 % NOTE: the behaviour of Fieldtrip's afni_niml_read is different - does not give cells 
 % can couse trouble, if Fieldtrips version is used
 if ~iscell(dannot_LH)
  dannot_LH{1}=dannot_LH;
 end
 if ~iscell(dannot_RH)
  dannot_RH{1}=dannot_RH;
 end

 if strcmp(dannot_LH{1}.nodes{1}.name,'SPARSE_DATA') && isa(dannot_LH{1}.nodes{1}.data,'double') % safety check
  annot_LH = dannot_LH{1}.nodes{1}.data; % assumes that annotation is in first node, usually the case
 else
  error('Unexpected trouble with reading the annotation file - LH')
 end
 if strcmp(dannot_LH{1}.nodes{3}.dset_type,'LabelTableObject') && isa(dannot_LH{1}.nodes{3}.nodes{1}.data,'cell') % safety check
  annot_key_LH = dannot_LH{1}.nodes{3}.nodes{1}.data(5:6); % assumes that annotation key is stored here (r,g,b,msk,index,label)
 else
  error('Unexpected trouble with reading the annotation key - LH')
 end

 if strcmp(dannot_RH{1}.nodes{1}.name,'SPARSE_DATA') && isa(dannot_RH{1}.nodes{1}.data,'double') % safety check
  annot_RH = dannot_RH{1}.nodes{1}.data; % assumes that annotation is in first node, usually the case
 else
  error('Unexpected trouble with reading the annotation file - RH')
 end
 if strcmp(dannot_RH{1}.nodes{3}.dset_type,'LabelTableObject') && isa(dannot_RH{1}.nodes{3}.nodes{1}.data,'cell') % safety check
  annot_key_RH = dannot_RH{1}.nodes{3}.nodes{1}.data(5:6); % assumes that annotation key is stored here (r,g,b,msk,index,label)
 else
  error('Unexpected trouble with reading the annotation key - RH')
 end

 % concatenate
 suma=[];
 suma.pos=vertcat(suma_LH.pos,suma_RH.pos);
 suma.tri=vertcat(suma_LH.tri,(suma_RH.tri+length(suma_LH.pos))); % shift trianage/faces count for RH
 suma.curv=vertcat(curv_LH,curv_RH);
 suma.annot=vertcat(annot_LH,annot_RH);
 suma.annot_key=cell(1,2);
 for i=1:2
  suma.annot_key{i}=vertcat(annot_key_LH{i},annot_key_RH{i});
 end 
 suma.msk=double(suma.annot~=0);
 suma.hemi=vertcat(repmat([1],length(suma_LH.pos),1),repmat([2],length(suma_LH.pos),1));
 suma.cortex=suma.msk; % here, all legitimate vertices are cortex

 % check length
 if (length(suma.curv) ~= length(suma.pos))
  error('mismatch between curvature and vertices, should not happen')
 end
 if (length(suma.annot) ~= length(suma.pos))
  error('mismatch between annotations and vertices, should not happen')
 end
end


end

