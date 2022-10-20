function [ subject ] = nmri_read_markers( subject, params )
%[ subject ] = nmri_read_markers( subject, params )
%   Will check spike/event markers and put into subject struct

% call the subject and params include
nmri_include_read_ps

%% Get the modality-specific analysis params
if isfield(subject,'dataset_mapping')
 [ params ] = nmri_get_dataset_params( params, subject.dataset_mapping );
end
[ params ] = nmri_get_modality_params( params, subject.dtype );

%% Safe subject in
subject_in=subject;

%% check the cleaned version
if (exist(subject.dws_filt_dataset,'file'))
 load(subject.dws_filt_dataset,'data');   
else
 error('We need the downsampled / filtered dataset to run marker reading. Use nmri_preproc first')
end

%% Mark spike trials (if this is called for by the analysis and we have events) - not actually reject them here
if (isfield(params,'rejectEvents') && params.rejectEvents == 1)
 
 % check for the event file
 if (isfield(params,'rejectEventsPath'))
  spath=params.rejectEventsPath;
 else
  spath=subject.analysis_dir;
 end
 if (isfield(params,'rejectEventsPattern'))
  pattern=strrep(params.rejectEventsPattern,'<subject_id>',subject.id);
 else
  pattern='*';
 end 
 if ischar(spath)
  spath={spath};
 end
 if ischar(pattern)
  pattern={pattern};
 end
 % check sizes
 if length(pattern)~=length(spath)
  error('Number of search paths and patterns mismatch - Check analysis_params')
 end
 % deal with raw_data paths
 spath=strrep(spath,'<raw_dataset>',subject.raw_dataset);
 spath=strrep(spath,'<analysis_dir>',subject.analysis_dir);
 % loop over paths
 found_files={};
 for iP=1:length(spath)
  eventFileThis=dir(fullfile(spath{iP},pattern{iP}));
  for iF=1:length(eventFileThis)
   found_files{end+1,1}=fullfile(eventFileThis(iF).folder,eventFileThis(iF).name);
  end
 end
 
 if (~isempty(found_files))
  % found at least one
  if (length(found_files)>1)
   % multiple, so choose
   manselect=listdlg('Name','Choose Event Marker File','PromptString',{'Found >1 possibility, please select (can be multiple)',['Subject=' subject.id  ', ExamID=' subject.exam_id ]},'SelectionMode','multiple','ListSize',[700 (50+(length(found_files)*10))],'ListString',found_files);
   if (~isempty(manselect))
    found_files=found_files(manselect);
   else
    error('No selection was made -- stopping here')
   end
  end
  
  if isfield(params,'rejectEventsTextscan')
   % Do textscan, usually fine for text files
   fid = fopen(found_files{1}); % file ID
   eval(params.rejectEventsTextscan);
   fclose(fid);
  elseif isfield(params,'rejectEventsFunction')
   % call a specific function, e.g. EGI xml oder EDF   
   events.evt_timings_seconds={};
   events.evt_IDs={};
   infos.evt_timings_seconds={};
   infos.evt_IDs={};
   for i=1:length(found_files)
    [ this_events, this_infos ]=feval(params.rejectEventsFunction,found_files{i},subject);
    events.evt_timings_seconds=[events.evt_timings_seconds; this_events.evt_timings_seconds];
    events.evt_IDs=[events.evt_IDs; this_events.evt_IDs];
    infos.evt_timings_seconds=[infos.evt_timings_seconds; this_infos.evt_timings_seconds];
    infos.evt_IDs=[infos.evt_IDs; this_infos.evt_IDs];
   end
   spikeTimes=cell2mat(events.evt_timings_seconds);
   spikeIDs=events.evt_IDs;
   infoTimes=cell2mat(infos.evt_timings_seconds);
   infoIDs=infos.evt_IDs;
  end
  
  spikeTimes = double(spikeTimes);
  nSpikes = length(spikeTimes);
  % check for spikeNames
  if ~exist('spikeIDs','var') || length(spikeIDs)~=nSpikes
   spikeIDs=repmat({'Spike'},nSpikes,1);
  end
  infoTimes = double(infoTimes);
  nInfos = length(infoTimes);
  if ~exist('infoIDs','var') || length(infoIDs)~=nInfos
   infoIDs=repmat({'Marking'},nInfos,1);
  end
                
  % convert the time into seconds of the spikes markers
  spikeTimesSeconds = spikeTimes/params.rejectEventsTimebase;
  infoTimesSeconds = infoTimes/params.rejectEventsTimebase;
                
  % sample number that has the spike markings
  spikes = round(spikeTimesSeconds*data.fsample)+1;
  infos = round(infoTimesSeconds*data.fsample)+1;
  
  subject.evt_markerFile=found_files{1};
  subject.evt_timings_sample=spikes;
  subject.evt_timings_seconds=spikeTimesSeconds;
  subject.evt_IDs=spikeIDs;
  
  subject.info_markerFile=found_files{1};
  subject.info_timings_sample=infos;
  subject.info_timings_seconds=infoTimesSeconds;
  subject.info_IDs=infoIDs;
  
  % redefine trials to discard the bad trials
  % the warning is okay
  % rd_cfg=[];
  % rd_cfg.trials = goodTrials;               
  % data = ft_redefinetrial(rd_cfg,data);
  % we do not remove them here, but consider this info later, NF 11/2016
  % 

  % check for (old) not found marker
  
  % if we got this far, place a stamp for completed event reading
  subject=nmri_stamp_subject(subject,'readingevents',params);
  if isfield(subject,'evt_markerFile_notFound')
   subject=rmfield(subject,'evt_markerFile_notFound');
  end
  
  fprintf('Found %d Events and %d Info Annotations\n',length(spikeIDs),length(infoIDs))
  
 else
  % no marker found, give warning if not a control
  if (~strcmp(subject.id(1:2),'C.'))
   h=msgbox({'No spike markings found: please make sure that this patient has no spikes!','',['Search path:' strjoin(spath,' ')],['Search pattern:' strjoin(pattern,' ')],'','Terminate processing, if needed'},'No events found');
   uiwait(h)
   subject.evt_markerFile_notFound=pattern;
  else
   disp('This seems to be a control subject and no marker file was found - Nothing was done')
  end
 end
else
 disp('No event rejection requested - Nothing was done')
end 

%% Parse stimuli if needed
dws_saved=0;
if isfield(params,'parse_stimuli') && params.parse_stimuli==1
 fprintf('Parseing stimuli...\n')
 data_out=nmri_parse_stimuli(subject, data, params);
 
 % and write out to dataset, if needed
 if ~isfield(data,'trial_markings') || ~isequal(data_out.trial_markings,data.trial_markings)
  fprintf('Saving new dws_filt dataset with stimuli info...\n')
  data=data_out;
  save(subject.dws_filt_dataset,'data','subject');
  dws_saved=1;
 end
end


%% Now write out things
if (~isequal(subject_in, subject))
 % save again, to be keep this info in mat file, unless done already
 if ~dws_saved
  save([subject.dws_filt_dataset],'subject','-append');
 end
 
 % add to other datasets, just in case
 if ~isfield(subject,'evt_markerFile_notFound') 
  % start with most advanced
  if isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')
   load(subject.hdm_lead,'subject')
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    % any change, so update
    subject.evt_markerFile=found_files{1};
    subject.evt_timings_sample=spikes;
    subject.evt_timings_seconds=spikeTimesSeconds;
    subject.evt_IDs=spikeIDs;
    subject.info_markerFile=found_files{1};
    subject.info_timings_sample=infos;
    subject.info_timings_seconds=infoTimesSeconds;
    subject.info_IDs=infoIDs;
    fprintf('Updating file = %s\n',subject.hdm_lead)
    save([subject.hdm_lead],'subject','-append');
   end
  end
  
  if isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')
   load(subject.cleanICA_dataset,'subject')
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    % any change, so update
    subject.evt_markerFile=found_files{1};
    subject.evt_timings_sample=spikes;
    subject.evt_timings_seconds=spikeTimesSeconds;
    subject.evt_IDs=spikeIDs;
    subject.info_markerFile=found_files{1};
    subject.info_timings_sample=infos;
    subject.info_timings_seconds=infoTimesSeconds;
    subject.info_IDs=infoIDs;
    fprintf('Updating file = %s\n',subject.cleanICA_dataset)
    save([subject.cleanICA_dataset],'subject','-append');
   end
  end
  
  if isfield(subject,'ICA_components') && exist(subject.ICA_components,'file')
   load(subject.ICA_components,'subject')    
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    % any change, so update
    subject.evt_markerFile=found_files{1};
    subject.evt_timings_sample=spikes;
    subject.evt_timings_seconds=spikeTimesSeconds;
    subject.evt_IDs=spikeIDs;
    subject.info_markerFile=found_files{1};
    subject.info_timings_sample=infos;
    subject.info_timings_seconds=infoTimesSeconds;
    subject.info_IDs=infoIDs;
    fprintf('Updating file = %s\n',subject.ICA_components)
    save([subject.ICA_components],'subject','-append');
   end
  end
  
  
  if isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
   load(subject.clean_dataset,'subject')    
   if ~isfield(subject,'evt_timings_seconds') || ~isequal(subject.evt_timings_seconds,spikeTimesSeconds) || ~isfield(subject,'info_timings_seconds') || ~isequal(subject.info_timings_seconds,infoTimesSeconds) ...
     || ~isfield(subject,'evt_IDs') || ~isequal(subject.evt_IDs,spikeIDs) || ~isfield(subject,'info_IDs')  || ~isequal( subject.info_IDs,infoIDs)
    % any change, so update
    subject.evt_markerFile=found_files{1};
    subject.evt_timings_sample=spikes;
    subject.evt_timings_seconds=spikeTimesSeconds;
    subject.evt_IDs=spikeIDs;
    subject.info_markerFile=found_files{1};
    subject.info_timings_sample=infos;
    subject.info_timings_seconds=infoTimesSeconds;
    subject.info_IDs=infoIDs;
    fprintf('Updating file = %s\n',subject.clean_dataset)
    save([subject.clean_dataset],'subject','-append');
   end
  end


 else
  % not file found -- still save
  if isfield(subject,'hdm_lead') && exist(subject.hdm_lead,'file')
   load(subject.hdm_lead,'subject')
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.hdm_lead)
    save([subject.hdm_lead],'subject','-append');
   end
  end 
  
  if isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file')
   load(subject.cleanICA_dataset,'subject')
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.cleanICA_dataset)
    save([subject.cleanICA_dataset],'subject','-append');
   end
  end 
  
  if isfield(subject,'ICA_components') && exist(subject.ICA_components,'file')
   load(subject.ICA_components,'subject')
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.ICA_components)
    save([subject.ICA_components],'subject','-append');
   end
  end 

  if isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file')
   load(subject.clean_dataset,'subject')    
   if ~isfield(subject,'evt_markerFile_notFound') || ~isequal(subject.evt_markerFile_notFound,pattern)
    % any change, so update
    subject.evt_markerFile_notFound=pattern;
    fprintf('Updating file = %s\n',subject.clean_dataset)
    save([subject.clean_dataset],'subject','-append');
   end
  end
 end
end
 


end

