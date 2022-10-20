function [ defaults ] = nf_request_spm( modality, spm_path )
%[ spm_defaults ] = nf_request_spm( modality )
%   Makes sure that SPM is in the path and defaults are sets
% modality and spm_path are optional (default, PET and ltest SPM12)

% firstoff, remove any SPM8's in the path
p=strsplit(path,pathsep);
for i=1:length(p)
 if (regexpi(p{i},'.*/spm8.*'))
  rmpath(p{i})
 end
end

if ~exist('modality','var')
 modality='PET';
end

if ~exist('spm_path','var')
 spm_path=fullfile(getenv('NMRI_TOOLS'),'/spm/spm12_7487'); % generic version
 if ~exist(spm_path,'dir')
  % fallback
  spm_path=fullfile(getenv('NMRI_TOOLS'),'/spm/spm12_7219'); %older version
  if ~exist(spm_path,'dir')
   % final fallback
   spm_path='/tools/spm/spm12_6685';
   if ~exist(spm_path,'dir')
    error('No SPM found, check paths')
   end
  end
 end
end
global defaults

if ~isfield(defaults,'modality')

 
 if isdeployed
  % either it is there, or we ware doomed.... ;)
  if exist('spm','file')==2
   fprintf('Setting up SPM in compiled mode\n')
   defaults = spm_get_defaults;
  else
   error('SPM was not included in the compiled archive, make sure to include in compile process. E.g. by setting cfg.spm=1 in nmri_qsub.m')
  end
 else
  if ~exist('spm','file')==2
   % usual mode, add to path if needed  
   addpath(spm_path);
  end
  spm('defaults',modality)
  fprintf('Setting up SPM in %s mode (non-compiled)\n',modality)
 end
 
 global defaults
 
end


end

