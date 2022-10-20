function [ subject_status, subject ] = nmri_check_subject_status( subject , params )
%[ subject_status, subject ] = nmri_check_subject_status( subject )
%   Will deterime if this subject is processed and return the flags

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
 
 % check the steps
 subject_status = [];
 % preproc
 if (isfield(subject,'dws_filt_dataset') && exist(subject.dws_filt_dataset,'file')) 
  subject_status.preproc=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'preproc'))
   subject_status.preproc_stamp=subject.stamps.preproc;
  end
 else
  subject_status.preproc=false;
 end
 

  
 % artifact rejection - visual
 if (isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')) 
  % donot consider ICA here
  subject_status.artifactrejection=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'artifactrejection'))
   subject_status.artifactrejection_stamp=subject.stamps.artifactrejection;
  end
 else
  subject_status.artifactrejection=false; 
 end
 
 
 
 %ICA estimation
 if (isfield(subject,'ICA_components') && exist(subject.ICA_components,'file')) 
  subject_status.artifactrejectionICA_estimate=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'artifactrejectionICA_estimate'))
   subject_status.artifactrejectionICA_estimate_stamp=subject.stamps.artifactrejectionICA_estimate;
  end
 else
  subject_status.artifactrejectionICA_estimate=false;
 end
 
  %ICA review
 if (isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')) 
  subject_status.artifactrejectionICA=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'artifactrejectionICA'))
   subject_status.artifactrejectionICA_stamp=subject.stamps.artifactrejectionICA;
  end
 else
  subject_status.artifactrejectionICA=false;
 end
 
 
 % Events
 if (isfield(subject,'evt_timings_seconds')) 
  subject_status.events=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'readingevents'))
   subject_status.events_stamp=subject.stamps.readingevents;
  end
 else
  % check if we want this
  if (isfield(params,'rejectEvents') && params.rejectEvents == 1)
   if ~strcmp(subject.id(1:2),'C.')
    if isfield(subject,'evt_markerFile_notFound')
     subject_status.events='none';
    else
     subject_status.events=false;
    end
   else
    subject_status.events='n. req';
   end
  else
   subject_status.events='n. req';
  end
 end
 
 % Vigilance
 if (isfield(params,'scoreVigilance') && params.scoreVigilance == 1)
  if isfield(subject,'data_info') && isfield(subject.data_info,'trial_markings')
   % we have some info, check if complete
   empty_matrix=cellfun(@(x) isempty(x),subject.data_info.trial_markings(:,:));
   % unscored items are empty in 1 and 4
   if size(empty_matrix,2)>3
    % we have provocations, these we do not score
    empty_matrix(~empty_matrix(:,4),1)=0;
   end
   % we also do not score bad trials
   if size(empty_matrix,2)>1
    % we have bad trials, these we do not score either
    empty_matrix(cellfun(@(x) islogical(x) && x==false,subject.data_info.trial_markings(:,2)),1)=0;
   end
   todo=sum(empty_matrix(:,1));
   if todo==0
    subject_status.vigilance=true;
   else
    subject_status.vigilance=sprintf('%d mis.',todo);
   end
   if (isfield(subject,'stamps') && isfield(subject.stamps,'vigilance'))
    subject_status.vigilance_stamp=subject.stamps.vigilance;
   end
  else
   subject_status.vigilance=false;
  end
 else
  subject_status.vigilance='n. req';
 end
 
 % MRI_ctf
 if (isfield(subject,'mri_ctf') && exist(subject.mri_ctf,'file')) 
  subject_status.mri_ctf=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'MRI_ctf'))
   subject_status.mri_ctf_stamp=subject.stamps.MRI_ctf;
  end
 else
  subject_status.mri_ctf=false; 
 end
 
 
 % MRI_seg
 if (isfield(subject,'mri_seg') && exist(subject.mri_seg,'file')) 
  subject_status.mri_seg=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'MRI_seg_SUMA'))
   subject_status.mri_seg_stamp=subject.stamps.MRI_seg_SUMA;
  end
 else
  subject_status.mri_seg=false; 
 end
 
 
 % Headmodel / Leadfield - need to know the headmodel for stamp
 if (isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')) 
  if (isfield(subject,'stamps') && isfield(params,'headmodel') && isfield(subject.stamps,['hdmSUMA_' params.headmodel]))
   subject_status.hdm_lead_stamp=subject.stamps.(['hdmSUMA_' params.headmodel]);
   subject_status.hdm_lead=true;
  else
   subject_status.hdm_lead=false;  
  end
 else
  subject_status.hdm_lead=false;   
 end
 
 % Processing - need to know the headmodel for stamp
 if (isfield(subject,'stats') && exist(subject.stats,'file')) 
  if (isfield(subject,'stamps') && isfield(params,'headmodel') && isfield(subject.stamps,['processing_' params.headmodel]))
   subject_status.processing_stamp=subject.stamps.(['processing_' params.headmodel]);
   subject_status.processing=true;
  else
   subject_status.processing=false; 
  end
 else
  subject_status.processing=false; 
 end
 
 % Processing Sensor
 if (isfield(subject,'sensor_stats') && exist(subject.sensor_stats,'file')) 
  subject_status.processing_sensor=true;
  if (isfield(subject,'stamps') && isfield(subject.stamps,'processing_sensor'))
   subject_status.processing_sensor_stamp=subject.stamps.('processing_sensor');
  end
 else
  subject_status.processing_sensor=false; 
 end
 
end

 
end

