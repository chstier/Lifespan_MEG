function nmri_check_seg(subject)
%nmri_check_seg(subject)
%   
%   Function to visually check (and edit if needed) the segmented MRI that
%   is the basis for the headmodel step.
%   Wraps around FSLeyes for the editing and deals with updating the
%   mri_seg.mat if needed

% check for fsleyes
[cstat,cout]=system('which fsleyes');

if cstat==0
 fsleyes_path=strtrim(cout);
else
 error('FSLeyes not in path. Make sure you add FSL to the system path before calling this script.')
end

if ~isfield(subject,'mri_seg') || ~exist(subject.mri_seg,'file')
 error('Cannot find mri_seg file in subject struct, or file is not available path. Run nmri_make_hdm_suma first.')
end


if ~isfield(subject,'mri_T1') || ~exist(subject.mri_T1,'file')
 error('Cannot find T1 image in subject struct, or file is not available path. Run nmri_align_mri first.')
end

load(subject.mri_seg)

[pa, fi, ~]=fileparts(subject.mri_seg);
work_path=fullfile(pa,'segmentation_visual');

if ~exist(work_path,'dir')
 mkdir(work_path)
end

% require SPM12 to be present
defaults=nf_request_spm;

% load the T1 as a space defining template
t1_vol=spm_vol(subject.mri_T1);

% now make tmp versions of the segmentations in mri_seg

labels={'tissue','skull','csf','white','grey'};

mT1=fullfile(work_path,['T1.nii']);
copyfile(subject.mri_T1,mT1);
setenv('FSLOUTPUTTYPE','NIFTI')
% get the transformation to std
[cstat, m2std]=system(['fslreorient2std -m ' fullfile(work_path,['m2std.mat']) ' ' mT1 ]);
% convert the matrix to numbers, need this to remap later
%m2std=strtrim(m2std);
%M=cellfun(@str2num,reshape(split(m2std),[4,4])');

% invert
system(['convert_xfm -omat ' fullfile(work_path,['m2orig.mat']) ' -inverse ' fullfile(work_path,['m2std.mat'])]);



% now re-orient for real
[cstat, m2std]=system(['fslreorient2std ' mT1 ' ' mT1]);

fslcall=[fsleyes_path ' ' mT1 ' -n T1 -cm greyscale' ];

% make a unified label image 
mri_u=zeros(mri_seg.dim);

for i=1:length(labels)
 mri_u(mri_seg.(labels{i}))=i;
end

% save unified
t1_vol.fname=fullfile(work_path,['all_classes.nii']);
if ~exist(t1_vol.fname,'file')
 fprintf('Writing unified NIFTI image\n')
 spm_write_vol(t1_vol,mri_u);
 system(['fslreorient2std ' t1_vol.fname ' ' t1_vol.fname]);
else
 fprintf('Unified NIFTI image present, delete if needed\n')
end
fslcall=[fslcall ' ' t1_vol.fname ' -ot label -a 50 --lut $NMRI_TOOLS/common/nmri-seg.lut'];

% now call FSLEYES

qcfg=[];
qcfg.question={'Now review the segmentation in FSLeyes',...
'Hit Alt-E to start editing and save the changes', 'Close/Quit FSLeyes when done.'};
qcfg.title='Review ';
qcfg.options={'Okay'};
qcfg.default={'Okay'};
qcfg.mode='buttons';
button=nf_uidialog(qcfg);

if strcmp(button,'Okay')
 system(fslcall);
 
 % check if you want to update
 qcfg=[];
 qcfg.question={'Do you want to update the segmentation with your changes?',...
 'Select Yes to overwrite old mri_seg or Abort/No to discard changes made now.'};
 qcfg.title='Keep edits?';
 qcfg.options={'Yes','No'};
 qcfg.default={'Yes'};
 qcfg.mode='buttons';
 button=nf_uidialog(qcfg);
 
 if strcmp(button,'Yes')
  % keep the edits and write back 
  modified=0;
  % transform back to native space
  system(['flirt -in ' fullfile(work_path,['all_classes.nii']) ' -out ' fullfile(work_path,['all_classes_back.nii']) ' -applyxfm -init ' fullfile(work_path,['m2orig.mat']) ' -ref ' subject.mri_T1]);
  mri_back=ft_read_mri(fullfile(work_path,['all_classes_back.nii']));
  
  for i=1:length(labels)
   fprintf('Checking class=%s\n',labels{i})
   tmp=mri_back.anatomy==i;
   if sum(mri_seg.(labels{i})~=tmp,'all')
    fprintf('N=%d voxels have been modified, updating class=%s\n',sum(mri_seg.(labels{i})~=tmp,'all'),labels{i})
    mri_seg.(labels{i})=tmp;
    modified=1;
   end
  end
  if modified==1
   fprintf('Saving new mri_seg\n\n')
   save(subject.mri_seg,'mri_seg');
   if isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')
    delete(subject.hdm_lead)
   else
    fprintf('Please make sure to delete any existing the hdm_lead file\n')
   end
   if isfield(subject,'electrodes_aligned') && exist(subject.electrodes_aligned,'file')
    delete(subject.electrodes_aligned)
   else
    fprintf('Please make sure to delete any existing the electrodes_aligned file (if EEG)\n')
   end
  else
   fprintf('No changes detected, keeping everything as-is\n\n')
  end
 end
 
 
 
end



end

