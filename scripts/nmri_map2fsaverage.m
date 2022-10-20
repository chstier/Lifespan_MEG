function [ all_metrics_fs_lh, all_metrics_fs_rh ] = nmri_map2fsaverage( all_subjects, all_metrics )
%[ all_metrics_fs_lh, all_metrics_fs_rh ] = nmri_map2fsaverage( all_subjects, all_metrics )
%   will transform any number of metric (cell array) to fsaverage

% Freesurfer defaults
suma_dir='/data/freesurfer/6.0.0/fsaverage/SUMA';
fs_dir='/data/freesurfer/6.0.0/fsaverage';

if nargin<2
 error('Input missing')
end

if isstruct(all_subjects)
 % in case a single subject is given
 all_subjects={all_subjects};
end

if ~iscell(all_subjects)
 error('Subject info needs to be cell array or struct (single case)')
end

if isnumeric(all_metrics)
 % in case a single subject is given
 all_metrics={all_metrics};
end

if ~iscell(all_metrics)
 error('Metrics needs to be cell array or numeric array (single case)')
end

% check N's
if length(all_metrics) ~= length(all_subjects)
 error('Mismatch of subjects and metrics number')
end

% get analyis_dir
analysis_dir=all_subjects{1}.analysis_dir;

N=length(all_subjects);

% get LD
SUMA_lds={};
SUMA_surfs={};
for n=1:N
 SUMA_lds{end+1}=all_subjects{n}.stamps.preproc.params.SUMA_ld;
 SUMA_surfs{end+1}=all_subjects{n}.stamps.preproc.params.SUMA_surf; 
end
SUMA_ld=unique(SUMA_lds);
SUMA_surf=unique(SUMA_surfs);

if length(SUMA_ld)~=1
 error('SUMA LD mismatch')
end
if length(SUMA_surf)~=1
 error('SUMA surf mismatch')
end
SUMA_ld=SUMA_ld{1};
SUMA_surf=SUMA_surf{1};

% now load the respective fsaverage - no subcortical here
fprintf('Loading Freesurfer-fsaverage surfaces\n')
file_LH = fullfile(fs_dir,'surf',['lh.' SUMA_surf]);
file_RH = fullfile(fs_dir,'surf',['rh.' SUMA_surf]);
if ~exist(file_LH,'file')
 error(['fsaverage surface not found = ' file_LH])
end
if ~exist(file_RH,'file')
 error(['fsaverage surface not found = ' file_RH])
end
[fs_LH.pos fs_LH.tri] = nf_load_surf_flexibel(file_LH);
[fs_RH.pos fs_RH.tri] = nf_load_surf_flexibel(file_RH);

% now load the SUMA template for our LD
load(fullfile(analysis_dir,'conf',['suma-all-fsaverage-' SUMA_ld '.mat']),'suma_all');

all_metrics_fs_lh=zeros(length(fs_LH.pos),1,1,N);
all_metrics_fs_rh=zeros(length(fs_RH.pos),1,1,N);
fprintf('Projecting from SUMA-ld %s to fsaverage...\n',SUMA_ld)
% save mapping for speed

% check files/dirs
mdir=fullfile(all_subjects{1}.analysis_dir,'conf','fsaverage-mappings');
if ~exist(mdir,'dir')
 mkdir(mdir);
end
mfile=fullfile(mdir,['mapping_SUMA_LD' SUMA_ld '_to' strrep(fs_dir,'/','_')  '.mat']);
if exist(mfile,'file')
 load(mfile,'nn_v_LH','nn_v_RH','nn_w_LH','nn_w_RH')
else
 nn_v_LH=cell(length(fs_LH.pos),1);
 nn_w_LH=cell(length(fs_LH.pos),1);
 nn_v_RH=cell(length(fs_RH.pos),1);
 nn_w_RH=cell(length(fs_RH.pos),1);
end
% now loop over all subjects
for n=1:N
 fprintf('%d/%d - subject: %s - LH',n,N,all_subjects{n}.id)
 
 if ~exist(mfile,'file') && sum(cellfun(@isempty,nn_v_LH))==0 && sum(cellfun(@isempty,nn_v_RH))==0
  % save if not there yet
  save(mfile,'nn_v_LH','nn_v_RH','nn_w_LH','nn_w_RH')
 end

 % now map the parameter - lh
 this_pos=(suma_all.hemi==1)&(suma_all.cortex<2); % choose hemi and no subcortex
 this_surf=[];
 this_surf.pos=suma_all.pos;
 this_surf.tri=suma_all.tri;
 % some need to go
 this_surf.pos(~this_pos,:)=repmat([NaN,NaN,NaN],sum(~this_pos),1);
 [this_surf.pos, this_surf.tri]=nf_remap_surface(this_surf.pos, this_surf.tri);
 % resample
 [all_metrics_fs_lh(:,1,1,n),~,nn_v_LH,nn_w_LH]=nf_surf2surf_transform(this_surf.pos,this_surf.tri,fs_LH.pos,all_metrics{n}(this_pos,1),nn_v_LH,nn_w_LH);
 
 fprintf(' - RH');
 % now map the parameter - rh
 this_pos=(suma_all.hemi==2)&(suma_all.cortex<2); % choose hemi and no subcortex
 this_surf=[];
 this_surf.pos=suma_all.pos;
 this_surf.tri=suma_all.tri;
 % some need to go
 this_surf.pos(~this_pos,:)=repmat([NaN,NaN,NaN],sum(~this_pos),1);
 [this_surf.pos, this_surf.tri]=nf_remap_surface(this_surf.pos, this_surf.tri);
 % resample
 [all_metrics_fs_rh(:,1,1,n),~,nn_v_RH,nn_w_RH]=nf_surf2surf_transform(this_surf.pos,this_surf.tri,fs_RH.pos,all_metrics{n}(this_pos,1),nn_v_RH,nn_w_RH);
 fprintf('...done\n')
end
fprintf('\n...all done\n')


end

