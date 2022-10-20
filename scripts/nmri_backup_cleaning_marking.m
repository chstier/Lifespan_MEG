function [ subject, backup ] = nmri_backup_cleaning_marking(subject, params)
%[ subject, backup ] = nmri_backup_cleaning_marking(subject, params)
%  
% This function will create a backup file of the cleaning, marking and vigilance
% scoring from one subject 
% Will use the dws-filt and cleaned/ICA-cleaned file to extract info
% Can be used to re-generate the cleaning
% Note: ICA-cleaning will have to be re-run 
% 

% subject       =   subject structure of cleaned subject, i.e. the one that was
%                   processed before.
%
% params        =   anaylsis parameter struct (optional, will search for
%                    analysis_params.m if not set)

% written by NF 10/2019


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct to work with')
end

% call the subject and params include
nmri_include_read_ps

% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


if (isfield(params,'useICA_clean') && params.useICA_clean==1)
 if (isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
  input=subject.cleanICA_dataset;
  useICA=true;
 else
  if (isfield(subject,'clean_dataset') &&  exist(subject.clean_dataset,'file'))
   input=subject.clean_dataset;
   useICA=false;
  else
   error('ICA cleaned and cleaned dataset not found - Run nmri_artifactrejection first')
  end
 end
else
 if (isfield(subject,'clean_dataset') &&  exist(subject.clean_dataset,'file'))
  input=subject.clean_dataset;
  useICA=false;
 else
  error('Cleaned dataset (w/o ICA) not found - Run nmri_artifactrejection first')
 end  
end
% retain subject info, if different
disp(['Loading cleaned dataset: ' input])
cleaned=load(input); 
if (~isfield(cleaned,'data') ) 
 error('Could not load cleaned data')
end

% safe some memory
cleaned.data.trial=[];


% now load the raw/unprocessed data
if (~isfield(subject,'dws_filt_dataset') || ~exist(subject.dws_filt_dataset,'file'))
 error('Filtered and downsampled dataset not specified, run nmri_preproc first')
else
 disp(['Loading raw dataset: ' subject.dws_filt_dataset ])
 raw=load(subject.dws_filt_dataset);
end
if (~isfield(raw,'data') ) 
 error('Could not load raw data')
end

% safe some memory
raw.data.trial=[];

% now extract

backup=[];
backup_date=datestr(now,'yyyymmdd');



% extract events - start  from cleaned
if isfield(cleaned.subject,'evt_timings_seconds')
 backup=backup_item(cleaned.subject,'evt_timings_sample',backup);
 backup=backup_item(cleaned.subject,'evt_timings_seconds',backup);
 backup=backup_item(cleaned.subject,'evt_IDs',backup);
elseif isfield(raw.subject,'evt_timings_seconds')
 backup=backup_item(raw.subject,'evt_timings_sample',backup);
 backup=backup_item(raw.subject,'evt_timings_seconds',backup);
 backup=backup_item(raw.subject,'evt_IDs',backup);
end


% make an empty trialmarkings
backup.trial_markings=cell(length(raw.data.time),4);
 % col1: sleep
 % col2: technical
 % col3: event
 % col4: rest/stimulation
 
backup.trial_markings_sampleinfo=cell(length(raw.data.time),2);
 % col1: sampleinfo start/stop
 % col2: seconds start/stop
for i=1:length(raw.data.time)
 backup.trial_markings_sampleinfo{i,1}=raw.data.sampleinfo(i,:);
 backup.trial_markings_sampleinfo{i,2}=raw.data.sampleinfo(i,:)/raw.data.fsample;
end



% compare cleaned and dws filt to fill trial_markings and bad channels for rejected trials 
for i=1:size(backup.trial_markings)
 % loop all original trials (there should never be a rejection here, even
 % with the oldest pipeline versions
 
 % match the time in cleaned
 maxmin=1; % want a 1 second match
 sel_trial=[];
 for ii=1:length(cleaned.data.time)
  minval=min(abs(cleaned.data.time{ii}-backup.trial_markings_sampleinfo{i,2}(1)));
  if minval<maxmin
   maxmin=minval;
   sel_trial=ii;
  end
 end
 
 if ~isempty(sel_trial)
  % check if we have a markup for this trial
  if isfield(cleaned.data,'trial_markings')
   % copy the marking
   backup.trial_markings(i,1:size(cleaned.data.trial_markings(sel_trial,:),2))=cleaned.data.trial_markings(sel_trial,:);
  else
   % no trial marking found, so nothing to do probably
  end
 else
  % could not find a match, so assume this trial was rejected/removed at
  % some point
  % this trials seems to have been rejected - mark as bad (false)
  backup.trial_markings{i,2}=false;  
 end
end

fprintf('Extracted valid trial markings for vigilance=%d\n',sum(~cellfun(@(x) isempty(x),backup.trial_markings(:,1))))
fprintf('Extracted valid trial markings for technical=%d\n',sum(~cellfun(@(x) isempty(x),backup.trial_markings(:,2))))
fprintf('Extracted valid trial markings for events=%d\n',sum(~cellfun(@(x) isempty(x),backup.trial_markings(:,3))))
fprintf('Extracted valid trial markings for stimulation=%d\n',sum(~cellfun(@(x) isempty(x),backup.trial_markings(:,4))))



% no deal with bad_channels
backup.bad_channels={};

if isfield(cleaned.data,'bad_channels')
 % these are bad for sure
 backup.bad_channels=cleaned.data.bad_channels;
end
% now look for rejected channels (pipeline version before 102019)
for i=1:length(raw.data.label)
 if ~any(strcmp(cleaned.data.label,raw.data.label{i}))
  % not found -seems bad
  if ~any(strcmp(backup.bad_channels,raw.data.label{i}))
   backup.bad_channels(end+1)=raw.data.label(i);
  end
 end
end
fprintf('Extracted bad_channnels, N=%d\n',length(backup.bad_channels))


% safe the backup

if ~isfield(subject,'backup')
 subject.backup=[];
end
if ~isfield(subject.backup,'date')
 % make new
 idx=1;
else
 % check idx
 idx=find(strcmp(subject.backup.date,backup_date));
 if length(idx)~=1
  % not found, create at the end
  idx=length(subject.backup)+1;
 end
end
subject.backup(idx).date=backup_date;
subject.backup(idx).data=backup;
subject.backup(idx).file=fullfile(subject.analysis_dir,subject.id,'backup',['backup_' subject.id '_' subject.exam_id '_' backup_date '.mat']);
if (~exist(fullfile(subject.analysis_dir,subject.id,'backup'),'dir'))
 mkdir(fullfile(subject.analysis_dir,subject.id,'backup'))
end

fprintf('Now saving the backup (%s)...\n',backup_date)
backup.date=backup_date;
save(subject.backup(idx).file,'backup');
fprintf('...done\n')



end

function backup=backup_item(orig,field,backup)
 if isfield(orig,field)
  backup.(field)=orig.(field);
 end
end


