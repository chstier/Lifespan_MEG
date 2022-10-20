function [ data ] = nmri_parse_stimuli( subject, data, params )
% will parse the info_ markings (from subject) to estimate trial-based
% stimuli markings
%
% Implemented for
% Augen auf / zu
% HV Beginn / End HV
% End HV / Post-HV 90 Sek
% ?? Hz / Aus


% call the subject and params include
nmri_include_read_ps

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% check the cleaned version
if ~exist('data','var') && exist(subject.dws_filt_dataset,'file')
 load(subject.dws_filt_dataset,'data');   
end
if ~exist('data','var')
 error('Either needs a data struct  or the downsampled / filtered dataset to run parsing. Use nmri_preproc first')
end

%% make the trial markings if needed
if ~isfield(data,'trial_markings') || size(data.trial_markings,1)==0
 data.trial_markings=cell(length(data.trial),4);
 % col1: sleep
 % col2: technical
 % col3: event
 % col4: rest/stimulation
 
 data.trial_markings_sampleinfo=cell(length(data.trial),2);
 % col1: sampleinfo start/stop
 % col2: seconds start/stop
 for i=1:length(data.trial)
  data.trial_markings_sampleinfo{i,1}=data.sampleinfo(i,:);
  data.trial_markings_sampleinfo{i,2}=data.sampleinfo(i,:)/data.fsample;
 end
end

if size(data.trial_markings,2)<4
 % make sure we have enough columns (backward compatability)
 for i=2:4
  if size(data.trial_markings,2)<i
   data.trial_markings=[data.trial_markings cell(size(data.trial_markings,1),1)];
  end
 end
end

if ~isfield(subject,'info_IDs') || ~isfield(subject,'info_timings_seconds')
 warning('No info markers in subject struct - cannot parse anything. Check if the markers have been read and re-run if needed')
else
%% Now do the actual parsing


 % PostHV
 startevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'End HV|HV Ende')));
 stopevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'Post HV 90')));
 grace_limit=460; % maximum time in seconds to accept
 grace_time=300; % if exceeding this time throw a warning and use grace_time instead
 data=parse_item(data,startevts,stopevts,grace_limit,grace_time,'postHV','Post-HV');
 
 % HV - will overwrite the postHV segment if needed
 startevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'HV Anfang|HV Beginn')));
 stopevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'End HV|HV Ende')));
 grace_limit=460; % maximum time in seconds to accept
 grace_time=300; % if exceeding this time throw a warning and use grace_time instead
 data=parse_item(data,startevts,stopevts,grace_limit,grace_time,'HV','HV');
 

 % eyes_open / closed - comes last, will overwrite HV
 startevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'Augen auf|eyes open')));
 stopevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'Augen zu|eyes closed')));
 grace_limit=60; % maximum time in seconds to accept
 grace_time=10; % if exceeding this time throw a warning and use grace_time instead
 data=parse_item(data,startevts,stopevts,grace_limit,grace_time,'eyes_open','Augen auf');
 
 
 % photostimulation will overwrite the rest
 startevts=subject.info_timings_seconds(~cellfun(@isempty,regexp(subject.info_IDs,'Hz$')));
 stopevts=subject.info_timings_seconds(~cellfun(@isempty,regexpi(subject.info_IDs,'^Off$')));
 grace_limit=60; % maximum time in seconds to accept
 grace_time=10; % if exceeding this time throw a warning and use grace_time instead
 data=parse_item(data,startevts,stopevts,grace_limit,grace_time,'PS','Photostimulation'); 
 
%%
end
%%

end

function data=parse_item(data,startevts,stopevts,grace_limit,grace_time,stim_ID,txt_label)
 
 if ~isempty(startevts)
  for ei=1:length(startevts)
   this_start=startevts(ei);
   this_stop=Inf;
   % find closest match after this
   rel_stop=stopevts-this_start;
   rel_stop(rel_stop<0)=[];
   rel_stop=sort(rel_stop);
  
   if ~isempty(rel_stop)
    % found at least one stop
    if (rel_stop(1))<grace_limit
     this_stop=rel_stop(1);
    else
     warning(sprintf('Could not find a legitimate stop marking for "%s" at %0.2fs. Next stop is at +%0.2fs, exceeding grace limit of %ds. Will use +%ds instead.',txt_label,this_start,rel_stop(1),grace_limit,grace_time))
     this_stop=grace_time;
    end
    this_stop=this_stop+this_start; % go back to absolute time again
   end
   
   % now find the right trial
   start_trial=[];
   stop_trial=[];
   for i=1:length(data.trial)
    if this_start>=data.time{i}(1,1) && this_start<=data.time{i}(1,end)
     start_trial=i;
    end
    if this_stop>=data.time{i}(1,1) && this_stop<=data.time{i}(1,end)
     stop_trial=i;
    end
   end
   
   if isempty(start_trial) && ~isempty(stop_trial)
    % start not in data, but stop - assume first trials
    start_trial=1;
   end
   if isempty(stop_trial) && ~isempty(start_trial)
    % stop not in data, but start - assume last trials
    stop_trial=length(data.trial);
   end
  
   if ~isempty(stop_trial) && ~isempty(start_trial) && start_trial<=stop_trial
    for i=start_trial:stop_trial
     data.trial_markings{i,4}=stim_ID;
    end
   end
  end
 else
  warning(sprintf('Could not dected any marking for "%s". This may be good and well, but could also indicate a problem. Investigate if needed.',txt_label))
 end
end
