function nf_batch_segment_new(files,outdir,def_field,bias,native,imported,unmodulated,modulated,cfg)
% function nf_batch_segment_new(files,outdir,def_field,bias,native,imported,unmodulated,modulated)
% Segements a batch of images with defaults otions
% Uses NEW segment (SPM8 + SPM12 only)
%
% files         = char array of files (if one channel)
%  ---or        = cell array of files (if multiple channels)
% outdir        = output dir 
% def_field     = 1: write y_ def field, 2: write iy_ def field, 3: write both y_ and iy_, 0: write nothing 
% native        = logical array of cX native images to write (e.g. [ 1 1 1 0 0 0] for c1,c2 and c3
% imported      = logical array of cX DARTEL imported images to write (e.g. [ 1 1 0 0 0 0] for rc1 and rc2
% unmodulated   = logical array of cX unmodulated, normalized images to write (e.g. [ 1 1 0 0 0 0] for wc1 and wc2
% modulated     = logical array of cX modulatedulated, normalized images to write (e.g. [ 1 1 0 0 0 0] for mwc1 and mwc2
% bias          = 1: write bias image, 2: write bias field, 3: write both, 0: write nothing 
% cfg           = struct for advanced configuration
% .bias_reg     = BIAS regularization (for SPM12 only), default: 0.001
% .bias_fwhm    = BIAS FMHM cutoff in mm(for SPM12 only), default: 60
%                 both BIAS params can be different per channel
%
% All auto-writes will be in 1.5mm resolution, for other write native and use deformation field
%
%
%
% written by Niels Focke 11/2010- 05/2017

global defaults

if ( ~isfield(defaults,'modality'))
 spm('PET');
end


spm('FnBanner','NF - Batch New Segment / Bias correct','2.20');

spm_progress_bar('Clear');


if ( ~exist('def_field','var') || ~isnumeric(def_field) || def_field<0 || def_field>4)
 def_field=spm_input('Write deformation fields ?','1','m','none...|forward...|inverse...|forward+inverse...',[0 1 2 3],2);
end

if ( ~exist('bias','var') || ~isnumeric(bias) || length(bias)~=1 || bias<0 || bias>3)
 bias=spm_input('Write BIAS corr image?','+1','m','none...|BIAS corr image...|BIAS corr field...|both...',[0 1 2 3],1);
end

while ( ~exist('native','var') || ~isnumeric(native) || length(native)~=6)
 native=spm_input('Write native c1...c6 (array)','+1','e',[1 1 1 1 1 0]);
end

while ( ~exist('imported','var') || ~isnumeric(imported) || length(imported)~=6)
 imported=spm_input('DARTEL imported rc1...rc6','+1','e',[1 1 0 0 0 0]);
end

while ( ~exist('unmodulated','var') || ~isnumeric(unmodulated) || length(unmodulated)~=6)
 unmodulated=spm_input('Normalized, unmodulated wc1...wc6','+1','e',[1 1 0 0 0 0]);
end

while ( ~exist('modulated','var') || ~isnumeric(modulated) || length(modulated)~=6)
 modulated=spm_input('Normalized, modulated mwc1...mwc6','+1','e',[1 1 0 0 0 0]);
end

if ~exist('cfg','var')
 cfg=[];
end
if ~isfield(cfg,'bias_reg') || ~isnumeric(cfg.bias_reg)
 cfg.bias_reg=0.001;
end
if ~isfield(cfg,'bias_fwhm') || ~isnumeric(cfg.bias_fwhm)
 cfg.bias_fwhm=60;
end

if ( ~exist('files','var') || isempty(files) )
 c=0;
 selectmore=1;
 while (selectmore==1)
  c=c+1;
  allfiles{c}=spm_select(inf, 'image', ['Choose images to be segmented (channel ' num2str(c) ') - ANY number']);
  if isempty(allfiles{c})
   error('Need to select at least one file');
  end
  selectmore=spm_input('Define an additional channel?','+1','m','no...|yes...',[0 1],0); 
 end
elseif (ischar(files))
 allfiles{1}=files;
elseif (iscell(files))
 allfiles=files;
end


% make BIAS params per channel
if length(allfiles)~=length(cfg.bias_reg)
 cfg.bias_reg=repmat(cfg.bias_reg(1),length(allfiles),1);
end
if length(allfiles)~=length(cfg.bias_fwhm)
 cfg.bias_fwhm=repmat(cfg.bias_fwhm(1),length(allfiles),1);
end

if ( ~exist('outdir','var') || isempty(outdir) )
 outdir=spm_select(1, 'dir', 'Choose destination dir (where output should be written to)');
end
if ( isempty(outdir) )
 outdir=[cd '/'];
end

% open log
fid=fopen([outdir '/segmentation_log'],'a+');
fprintf(fid,'Start Time:\t%s\n',datestr(now));

%-----------------------------------------------------------------------
if (exist('defaults.preproc'))
 % spm8 tree
 for i=1:length(allfiles)
  cellfiles{i}=cellstr(allfiles{i});
  matlabbatch{1}.spm.tools.preproc8.channel(i).vols = cellfiles{i};
  matlabbatch{1}.spm.tools.preproc8.channel(i).biasreg = cfg.bias_reg(i);
  matlabbatch{1}.spm.tools.preproc8.channel(i).biasfwhm = cfg.bias_fwhm(i); 
 end
 
 matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = { [ spm('dir') '/toolbox/Seg/TPM.nii,1' ] };
 matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
 
 matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = { [ spm('dir') '/toolbox/Seg/TPM.nii,2' ] };
 matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
 
 matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = { [ spm('dir') '/toolbox/Seg/TPM.nii,3' ] };
 matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
 
 matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = { [ spm('dir') '/toolbox/Seg/TPM.nii,4' ] };
 matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3; 
 
 matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = { [ spm('dir') '/toolbox/Seg/TPM.nii,5' ] };
 matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4; 
 
 matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = { [ spm('dir') '/toolbox/Seg/TPM.nii,6' ] };
 matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
 
 matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
 matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
 matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;

 for i=1:6
  matlabbatch{1}.spm.tools.preproc8.tissue(i).native = [native(i) imported(i)]; 
  matlabbatch{1}.spm.tools.preproc8.tissue(i).warped = [unmodulated(i) modulated(i)];
 end
 
 for i=1:length(allfiles)
  switch bias
   case 0
   matlabbatch{1}.spm.tools.preproc8.channel(i).write = [0 0]; % no bias
   case 1
   matlabbatch{1}.spm.tools.preproc8.channel(i).write = [0 1]; % bias image
   case 2
   matlabbatch{1}.spm.tools.preproc8.channel(i).write = [1 0]; % bias field
   case 3
   matlabbatch{1}.spm.tools.preproc8.channel(i).write = [1 1]; % both
  end
 end
  
 switch def_field
  case 0
  matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0]; % no def field
  case 1
  matlabbatch{1}.spm.tools.preproc8.warp.write = [0 1]; % forward def field
  case 2
  matlabbatch{1}.spm.tools.preproc8.warp.write = [1 0]; % inverse def field
  case 3
  matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1]; % both def fields
 end
 
 fprintf('BIAS regularization (channel %d):\t%f\n',i,matlabbatch{1}.spm.tools.preproc8.channel.biasreg);
 fprintf('BIAS FWHM (channel %d):\t%f\n',i,matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm);
 
 % and print to file 
 fprintf(fid,'BIAS FWHM (channel %d):\t%f\n',i,matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm);
 fprintf(fid,'BIAS regularization (channel %d):\t%f\n',i,matlabbatch{1}.spm.tools.preproc8.channel.biasreg);

else
 % SPM 12 tree
 for i=1:length(allfiles)
  cellfiles{i}=cellstr(allfiles{i});
  matlabbatch{1}.spm.spatial.preproc.channel(i).vols = cellfiles{i};
  matlabbatch{1}.spm.spatial.preproc.channel(i).biasreg = cfg.bias_reg(i);
  matlabbatch{1}.spm.spatial.preproc.channel(i).biasfwhm = cfg.bias_fwhm(i);
  fprintf('BIAS regularization (channel %d):\t%f\n',i,matlabbatch{1}.spm.spatial.preproc.channel(i).biasreg);
  fprintf('BIAS FWHM (channel %d):\t%f\n',i,matlabbatch{1}.spm.spatial.preproc.channel(i).biasfwhm);
  % and print to file 
  fprintf(fid,'BIAS regularization (channel %d):\t%f\n',i,matlabbatch{1}.spm.spatial.preproc.channel(i).biasreg);
  fprintf(fid,'BIAS FWHM (channel %d):\t%f\n',i,matlabbatch{1}.spm.spatial.preproc.channel(i).biasfwhm);
 end
 
 matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = { [ spm('dir') '/tpm/TPM.nii,1' ] };
 matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1; %this is also different in SPM12
 
 matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = { [ spm('dir') '/tpm/TPM.nii,2' ] };
 matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1; %this is also different in SPM12
 
 matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = { [ spm('dir') '/tpm/TPM.nii,3'  ] };
 matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
 
 matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = { [ spm('dir') '/tpm/TPM.nii,4'  ] };
 matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3; 
 
 matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = { [ spm('dir') '/tpm/TPM.nii,5' ] };
 matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4; 
 
 matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = { [ spm('dir') '/tpm/TPM.nii,6' ] };
 matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
 
 matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1; % SPM12 default, was of/not present in SPM8
 matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
 matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];

 matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
 matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;

 for i=1:6
  matlabbatch{1}.spm.spatial.preproc.tissue(i).native = [native(i) imported(i)]; 
  matlabbatch{1}.spm.spatial.preproc.tissue(i).warped = [unmodulated(i) modulated(i)];
 end
 
 for i=1:length(allfiles)
  switch bias
   case 0
   matlabbatch{1}.spm.spatial.preproc.channel(i).write = [0 0]; % no bias
   case 1
   matlabbatch{1}.spm.spatial.preproc.channel(i).write = [0 1]; % bias image
   case 2
   matlabbatch{1}.spm.spatial.preproc.channel(i).write = [1 0]; % bias field
   case 3
   matlabbatch{1}.spm.spatial.preproc.channel(i).write = [1 1]; % both
  end
 end
  
 switch def_field
  case 0
  matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0]; % no def field
  case 1
  matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1]; % forward def field
  case 2
  matlabbatch{1}.spm.spatial.preproc.warp.write = [1 0]; % inverse def field
  case 3
  matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % both def fields
 end
 
 
end
save([outdir '/segmentation.mat'],'matlabbatch');

spm_jobman('run',matlabbatch);


for i=1:length(cellfiles{1})
 [pa, fi, ext] = fileparts(cellfiles{1}{i});
 fprintf ([ 'Image = ' fi ' \n'])
 fprintf(fid,[ 'Image = ' fi ' \n']);

 if ( ~strcmp(pa,outdir) ) && ( ~strcmp([pa '/'],outdir) ) 
  % move classes
  for c=1:6
   if (native(c)==1)
    movefile([ pa '/c' num2str(c) fi '*' ],outdir);
   end
   if (imported(c)==1)
    movefile([ pa '/rc' num2str(c) fi '*' ],outdir);
   end
   if (unmodulated(c)==1)
    movefile([ pa '/wc' num2str(c) fi '*' ],outdir);
   end
   if (modulated(c)==1)
    movefile([ pa '/mwc' num2str(c) fi '*' ],outdir);
   end  
  end
  
  movefile([ pa '/' fi '_seg8.mat' ],outdir);
  
  oldpa=pa;
  
  for ii=1:length(cellfiles)
   [pa, fi, ext] = fileparts(cellfiles{ii}{i});
   % BIAS images
   switch bias
   case 1
    movefile([ pa '/m' fi '*' ],outdir);
   case 2
    movefile([ pa '/BiasField_' fi '*' ],outdir);
   case 3
    movefile([ pa '/m' fi '*' ],outdir);
    movefile([ pa '/BiasField_' fi '*' ],outdir);
   end
  end
   
  % Def images - handels differently for SPM8 and SPM12
  if (exist('defaults.preproc'))
   % SPM8 tree, def field is named after the last class
   [pa, fi, ext] = fileparts(cellfiles{end}{i});
  else
   %SPM12 takes first cell   
   [pa, fi, ext] = fileparts(cellfiles{1}{i});
  end
  
  switch def_field
  case 1
   movefile([ oldpa '/y_' fi '*' ],outdir);
  case 2
   movefile([ oldpa '/iy_' fi '*' ],outdir);
  case 3
   movefile([ oldpa '/y_' fi '*' ],outdir);
   movefile([ oldpa '/iy_' fi '*' ],outdir);
  end
  
 end
end
fclose(fid);

