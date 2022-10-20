function [ all_subjects ] = nmri_all_subjects( analysis_dir, use_cache )
%[ all_subjects ] = nmri_all_subjects( analysis_dir, use_cache)
%  Load all subject structs that are in the current analyis dir, will also
%  re-root moved structs

if ~exist('use_cache','var') || isempty(use_cache)
 use_cache=0; % default to always read from disk
end

% check the params struct
if (~exist('params','var'))
 prev_dir=pwd;
 cd(analysis_dir)
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in analysis dir ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the paramter file (analysis_params.m)')  
  end
 end
 cd(prev_dir)
end

% set cache dir, if requested
cache=[];
cache.id={};
cache.subjects={};
cache.timestamp=[];
if use_cache
 cache_file=fullfile(analysis_dir,'subjects_cache.mat');
 if exist(cache_file,'file')
  load(cache_file,'cache')
 end
end


all_dirs=dir(analysis_dir);
all_subjects={};

for i=1:size(all_dirs,1)
 % check for conf dir subject_info.m or subject_info.mat
 if ~strcmp(all_dirs(i).name(1),'.') && all_dirs(i).isdir && exist(fullfile(analysis_dir,all_dirs(i).name,'conf'),'dir')
  % have conf dir, now look for subject infos
  all_exams=dir(fullfile(analysis_dir,all_dirs(i).name,'conf'));
  %tic
  %disp(all_dirs(i).name)
  for ii=1:size(all_exams,1)
   if ~strcmp(all_exams(ii).name(1),'.')
    % use .mat file on most cases
    if exist(fullfile(analysis_dir,all_dirs(i).name,'conf',all_exams(ii).name,'subject_info.mat'),'file')
     % check if this is newer than cache
     mfile=dir(fullfile(analysis_dir,all_dirs(i).name,'conf',all_exams(ii).name,'subject_info.mat'));
     cid=find(strcmpi(cache.id,[all_dirs(i).name '-' all_exams(ii).name]));
     if isempty(cid) || cache.timestamp(cid)<mfile.datenum 
      % not cached or newer
      load(fullfile(analysis_dir,all_dirs(i).name,'conf',all_exams(ii).name,'subject_info.mat'),'subject');
     elseif ~isempty(cid)
      % use cached version
      subject=cache.subjects{cid};
     end 
     
     if ~isfield(subject,'exam_id')
      subject.exam_id=all_exams(ii).name;
     end
     if (~strcmp(subject.analysis_dir,analysis_dir))
      % check now for move and remember to check again later
      subject=nmri_reroot_subject(subject,analysis_dir); 
      reroot=1;
     else
      reroot=0;
     end
     
     % now check if we have a cached version
     do_load=1;
     if use_cache 
      cid=find(strcmpi(cache.id,[subject.id '-' subject.exam_id]));
      if ~isempty(cid) && isfield(subject,'loaded_from') && exist(subject.loaded_from,'file')
       mfile=dir(subject.loaded_from);
       if cache.timestamp(cid)>=mfile.datenum
        % we have a potential cache hit
        
        % check for newer files downstream
        if (~isfield(subject,'SelectedTrials_file'))
         subject.SelectedTrials_file=fullfile(subject.analysis_dir,subject.id,'processed',['selected_trials_' subject.id '_' subject.exam_id '.mat']);
        end
        if (~isfield(subject,'hdm_lead') && isfield(params,'headmodel'))
         subject.hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
        end
        if (~isfield(subject,'cleanICA_dataset'))
         subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
        end
        if (~isfield(subject,'clean_dataset'))
         subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
        end
        if (~isfield(subject,'dws_filt_dataset'))
         subject.dws_filt_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_' subject.id '_' subject.exam_id '.mat']);
        end
        
        check_all={subject.SelectedTrials_file,subject.cleanICA_dataset,subject.clean_dataset,subject.dws_filt_dataset};
        if isfield(params,'headmodel')
         check_all=[check_all subject.hdm_lead];
        end
        
        do_load=0;
        for ci=1:length(check_all)
         if exist(check_all{ci},'file')
          mfile=dir(check_all{ci});
          if cache.timestamp(cid)<mfile.datenum
           % newer file, so load
           do_load=1;
          end
         end
        end
        
        % we still have a valid hit, so read from cache
        if ~do_load
         this_subject=cache.subjects{cid};
        end
        
       end
      end
     end
     % now load the most advanced, if not cached / hit
     if do_load
      fprintf('Loading from disk for: %s-%s\n',subject.id,subject.exam_id)
      this_subject=nmri_load_subject_most_advanced(subject,params);
      % and save to cache
      if use_cache
       cid=find(strcmpi(cache.id,[this_subject.id '-' this_subject.exam_id]));
       if isempty(cid)
        cache.id{end+1}=[this_subject.id '-' this_subject.exam_id];
        cache.subjects{end+1}=this_subject;
        mfile=dir(this_subject.loaded_from);
        cache.timestamp(end+1)=mfile.datenum;
       else
        cache.subjects{cid}=this_subject;
        mfile=dir(this_subject.loaded_from);
        cache.timestamp(cid)=mfile.datenum;
       end
      end
     end
     % finally, check if this subject struct was moved, i.e. if there is a 
     % mismatch of analyis_dirs
     if reroot
      this_subject.rerooted=sprintf('Original Dir: %s, rerooted to: %s',this_subject.analysis_dir,analysis_dir);
      this_subject=nmri_reroot_subject(this_subject,analysis_dir); 
     end
     all_subjects{end+1,1}=this_subject;
    
    % legacy subject_info.mat
    elseif exist(fullfile(analysis_dir,all_dirs(i).name,'conf',all_exams(ii).name,'subject_info.m'),'file') 
     if ~isdeployed
      cd(fullfile(analysis_dir,all_dirs(i).name,'conf',all_exams(ii).name))
      clear subject
      eval('subject_info');
      cd(analysis_dir);
     else
      error('Cannot use eval in deployed/compiled mode')
     end
     
     if ~isfield(subject,'exam_id')
     subject.exam_id=all_exams(ii).name;
     end
     if (~strcmp(subject.analysis_dir,analysis_dir))
      % check now for move and remember to check again later
      subject=nmri_reroot_subject(subject,analysis_dir); 
      reroot=1;
     else
      reroot=0;
     end
     this_subject=nmri_load_subject_most_advanced(subject,params);
     % finally, check if this subject struct was moved, i.e. if there is a 
     % mismatch of analyis_dirs
     if reroot
      this_subject.rerooted=sprintf('Original Dir: %s, rerooted to: %s',this_subject.analysis_dir,analysis_dir);
      this_subject=nmri_reroot_subject(this_subject,analysis_dir); 
     end
     all_subjects{end+1,1}=this_subject;
    end
    
   end
  end
  %toc
 end

end



% now check for group (dx_filter) files - can be done 
% in analyis dir
all_dirs=dir(fullfile(analysis_dir,'*dx_filter*'));
dx_files={};
for i=1:size(all_dirs,1)
 if all_dirs(i).bytes>0
  dx_files{end+1}=fullfile(analysis_dir,all_dirs(i).name);
 end
end

% maybe also in group_filters
if exist(fullfile(analysis_dir,'group_filters'),'dir')
 all_dirs=dir(fullfile(analysis_dir,'group_filters','*dx_filter*'));
 for i=1:size(all_dirs,1)
  if all_dirs(i).bytes>0
   dx_files{end+1}=fullfile(analysis_dir,'group_filters',all_dirs(i).name);
  end
 end
end

if ~isempty(dx_files)
 % found at least one - read
 all_cols={};
 for i=1:length(dx_files)
  fid=fopen(dx_files{i},'r');
  cols=textscan(fid,'%s %s %*[^\n]','HeaderLines',1);
  all_cols=vertcat(all_cols,cols);
  fclose(fid);
 end
 % now parse all subject
 for i=1:length(all_subjects)
  % check by regular expression
  fdx=~cellfun(@isempty,regexpi(all_cols{1},[ '^' all_subjects{i}.id '.*' ]));
  if sum(fdx)>0
   grp=unique(all_cols{2}(fdx));
   if length(grp)>1
    warning(['Found >1 group match for ID=' all_subjects{i}.id '. Will take only the first one( ' grp{1} ' ). Still, check your group files...'])
   end   
   all_subjects{i}.group=grp{1};
  else
   % nothing found, make one based in ID/basname
   if strcmp(all_subjects{i}.id(1),'C')
    all_subjects{i}.group='Cx';
   elseif  strcmp(all_subjects{i}.id(1),'P')
    all_subjects{i}.group='Pat';
   end
  end
 end
else
 fprintf('No group filters (dx_groups) found. No group field will be set.\n') 
end

% save the cache file, if requested
if use_cache
 save(cache_file,'cache')
end

end

