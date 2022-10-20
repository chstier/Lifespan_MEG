function [ subject ] = nmri_load_subject_most_advanced( subject, params )
%[ subject ] = nmri_load_subject_most_advanced( subject )
%   Will load the most advanced subject struct based on step-hirarchy

% check if subject struct is okay
if ~isfield(subject,'id') || ~isfield(subject,'exam_id') || ~isfield(subject,'analysis_dir') || ~isfield(subject,'dtype')
 subject_status = false;
 warning(['Subject struct is lacking basic info -- should not happen'])
else
 % this seems to be at least basically correct
 
 % check the params struct
 if (~exist('params','var'))
  prev_dir=pwd;
  cd(subject.analysis_dir)
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
 [ params ] = nmri_get_modality_params( params, subject.dtype );
 
 % preserve pre-load analyis dir
 analysis_dir=subject.analysis_dir;

 % load from selected trials, if there
 if (~isfield(subject,'SelectedTrials_file'))
  subject.SelectedTrials_file=fullfile(subject.analysis_dir,subject.id,'processed',['selected_trials_' subject.id '_' subject.exam_id '.mat']);
 end
 if (exist(subject.SelectedTrials_file,'file'))
  load_file=subject.SelectedTrials_file;
  load(load_file,'subject');
  subject.loaded_from=load_file;
 else
  % not there, check for hdm_lead
  if (~isfield(subject,'hdm_lead') && isfield(params,'headmodel'))
   subject.hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
  end
  if (isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file'))
   load_file=subject.hdm_lead;
   load(load_file,'subject');
   subject.loaded_from=load_file;
  else
   % not there, check for ICA cleaned datset
   if (~isfield(subject,'cleanICA_dataset'))
    subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
   end
   if (exist(subject.cleanICA_dataset,'file'))
    load_file=subject.cleanICA_dataset;
    load(load_file,'subject');
    subject.loaded_from=load_file;   
   else 

    % not there, check for artifact_corrected
    if (~isfield(subject,'clean_dataset'))
     subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
    end
    if (exist(subject.clean_dataset,'file'))
     load_file=subject.clean_dataset;
     load(load_file,'subject');
     subject.loaded_from=load_file; 
    else
     % not there, then check fo dws_filt
     if (~isfield(subject,'dws_filt_dataset'))
      subject.dws_filt_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_' subject.id '_' subject.exam_id '.mat']);
     end
     if (exist(subject.dws_filt_dataset,'file'))
      load_file=subject.dws_filt_dataset;
      load(load_file,'subject');
      subject.loaded_from=load_file; 
     else
      % not even this is present - load the .m file
      mdir=fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id);
      if exist(mdir,'file')
       if exist(fullfile(mdir,'subject_info.mat'),'file')
        load(fullfile(mdir,'subject_info.mat'),'subject');
        subject.loaded_from=fullfile(mdir,'subject_info.mat');
       elseif exist(fullfile(mdir,'subject_info.m'),'file')
        prev_dir=pwd;
        cd(mdir)
        eval('subject_info');
        cd(prev_dir)
        subject.loaded_from=fullfile(mdir,'subject_info.m');
       end
      end
     end
    end
   end
  end
 end

 % deal with reroot
 if (~strcmp(subject.analysis_dir,analysis_dir))
  subject=nmri_reroot_subject(subject,analysis_dir); 
 end
 % get rid of data sometimes stored in the hdr for MFF
 if isfield(subject,'hdr') && isfield(subject.hdr,'orig') && isfield(subject.hdr.orig,'data')
  subject.hdr.orig=rmfield(subject.hdr.orig,'data');
 end
  if isfield(subject,'hdr') && isfield(subject.hdr,'orig') && isfield(subject.hdr.orig,'event')
  subject.hdr.orig=rmfield(subject.hdr.orig,'event');
 end
  
 % special jiggle for MRI_ctf, in case it is processed but stamp not loaded
 if ((~isfield(subject,'stamps') || ~isfield(subject.stamps,'MRI_ctf')) && exist(fullfile(subject.analysis_dir,subject.id,'logs',subject.exam_id,'MRI_ctf.mat'),'file'))
  t=load(fullfile(subject.analysis_dir,subject.id,'logs',subject.exam_id,'MRI_ctf.mat'),'subject');
  subject.stamps.MRI_ctf=t.subject.stamps.MRI_ctf;
 end
 % make sure we always have some relevant fields
 if (~isfield(subject,'mri_dir'))
  subject.mri_dir=fullfile(subject.analysis_dir,subject.id,'mri');
 end
 if (~isfield(subject,'mri_ctf'))
  subject.mri_ctf=fullfile(subject.mri_dir,['mri_ctf_' subject.id '.mat']);
 end
 if (~isfield(subject,'mri_seg'))
  subject.mri_seg=fullfile(subject.mri_dir,['mri_seg_' subject.id '.mat']);
 end
 if (~isfield(subject,'suma_surface'))
  subject.suma_surface=fullfile(subject.mri_dir,['suma_surface_' subject.id '.mat']);
 end
 if (~isfield(subject,'clean_dataset'))
  subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
 end
 if (isfield(params,'useICA_clean') && params.useICA_clean==1)
  if (~isfield(subject,'cleanICA_dataset'))
   subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
  end
 end
 if (~isfield(subject,'stats_dir'))
  subject.stats_dir=fullfile(subject.analysis_dir,subject.id,'stats');
 end
 if (~isfield(subject,'sensor_stats') && isfield(params,'headmodel'))
  subject.sensor_stats=fullfile(subject.stats_dir,['sensor_stats_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
 end
 
 % deal with missing stamps
 if ~isfield(subject,'stamp_dir')
  subject.stamp_dir=fullfile(subject.analysis_dir,subject.id,'logs',subject.exam_id);
 end
 if (~isfield(subject,'stamps'))
  subject.stamps=[];
 end
 items={'preproc','artifactrejection','artifactrejectionICA','MRI_ctf','MRI_seg_SUMA','vigilance','processing_sensor'};
 if isfield(params,'headmodel')
  items=[items {['hdmSUMA_' params.headmodel],['processing_' params.headmodel]}];
 end
 for i=1:length(items)
  if ~isfield(subject.stamps,items{i}) && exist(fullfile(subject.stamp_dir,[items{i} '.mat']),'file')
   % load missed stamp
   this_stamp=load(fullfile(subject.stamp_dir,[items{i} '.mat']),'subject');
   if isfield(this_stamp.subject.stamps,items{i})
    subject.stamps.(items{i})=this_stamp.subject.stamps.(items{i});
   %for some time we had a different place for stamps...
   elseif isfield(this_stamp.subject.stamps,subject.exam_id) && isfield(this_stamp.subject.stamps.(subject.exam_id),items{i})
    subject.stamps.(items{i})=this_stamp.subject.stamps.(subject.exam_id).(items{i});
   end
  end   
 end
 
end

end

