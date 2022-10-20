function [ goodTrials, badTrials ] = nmri_trial_selector( subject, data, params, output )
% [ goodTrials, badTrials ] = nmri_trial_selector( subject, data, params )
%   Will determine god and bad trials, based on the selectors in
%   analysis_params, trial_info and events
%   will automatically be called during trial selection step
%   (nmri_processing/_sensor)



% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

if (~exist('output','var') ) 
 output=true;
end

% call the subject and params include
nmri_include_read_ps


if ~exist('data','var') || isempty(data) 
 % check if we have cleaned dataset (ICA or not)
 if (isfield(params,'useICA_clean') && params.useICA_clean==1)
  if (isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
   input=subject.cleanICA_dataset;
   useICA=true;
  else
   error('ICA cleaned dataset not found - Run nmri_artifactrejection first')
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
 load(input,'data'); 
 if (~exist('data','var') ) 
  error('Could not load data')
 end
end

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );
% and the mapping, if availabel
if isfield(subject,'dataset_mapping')
 [ params ] = nmri_get_dataset_params( params, subject.dataset_mapping );
end


%% No to the selection
if (isfield(params,'nTrials'))
 nTrials=params.nTrials;
else
 nTrials=[]; % take all that are good
end

% determine trials now
badTrials = [];
if isfield(data,'ntrials') % true for data_info
 dataNtrials=data.ntrials;
else
 % start with all good
 dataNtrials=length(data.trial);
end
goodTrials = [1:dataNtrials];


 
if (isfield(params,'rejectEvents') && params.rejectEvents == 1)
 % check in dws_filt dataset not to miss atypically run processing
 if ~isfield(subject,'evt_timings_seconds') && ~isfield(subject,'evt_markerFile_notFound')
  subject_clean=load(subject.dws_filt_dataset,'subject');
  items=fieldnames(subject_clean.subject);
  for i=1:length(items)
   if length(items{i})>4 && strcmp(items{i}(1:4),'evt_')
    subject.(items{i})=subject_clean.subject.(items{i});
   end
  end
  if isfield(subject_clean.subject,'stamps') && isfield(subject_clean.subject.stamps,'readingevents')
   subject.stamps.readingevents=subject_clean.subject.stamps.readingevents;
  end
  clear subject_clean
 end

 if isfield(subject,'evt_timings_seconds')
  % we have events - so reject these trials
  nSpikes=length(subject.evt_timings_seconds);
  trial_samples=round(params.trial_length*data.fsample);
  
  % check for events
  for iTrial = 1:dataNtrials
   for iSpike = 1:nSpikes
    if (data.time{iTrial}(1,end) > subject.evt_timings_seconds(iSpike)) && (data.time{iTrial}(1,1) < subject.evt_timings_seconds(iSpike))
     if output
      disp(['bad trial (by event): ' num2str(iTrial)]);
     end
     badTrials = [badTrials iTrial];
     for remTi=1:length(params.rejectEventsTrial)
      % now look at the neighbours
      remTrial=iTrial+params.rejectEventsTrial(remTi);
      if (remTrial~=iTrial) && (remTrial > 0) && (remTrial <= dataNtrials) && ...
        (data.time{remTrial}(1,end) > subject.evt_timings_seconds(iSpike) +(params.rejectEventsTrial(remTi)*params.trial_length)) && ...
        (data.time{remTrial}(1,1) < subject.evt_timings_seconds(iSpike)+(params.rejectEventsTrial(remTi)*params.trial_length))
       badTrials = [badTrials remTrial];
       if output
        disp(['bad trial (by event neighbour): ' num2str(remTrial)]);
       end
      end
     end
    end
   end
  end
 
 end
end % done with event-based rejection

% now check for trial markings 
if isfield(data,'trial_markings')
 remTrial=find(cellfun(@(x) islogical(x) && x==false,data.trial_markings(:,2)))';
 if output
  fprintf('Found %d bad trial(s) based on technical marking\n',length(remTrial))
 end
 badTrials = [badTrials remTrial];
end
if isfield(params,'rejectVigilance') && isstruct(params.rejectVigilance) && isfield(data,'trial_markings')
 % reject based on vigilance
 items=fieldnames(params.rejectVigilance);
 for i=1:length(items)
  if params.rejectVigilance.(items{i})==1
   % reject this one
   remTrial=find(strcmp(data.trial_markings(:,1),items{i}))';
   if output
    fprintf('Found %d bad trial(s) based on vigilance=%s\n',length(remTrial),items{i})
   end
   badTrials = [badTrials remTrial];
  end
 end
end
if isfield(params,'rejectStimuli') && isstruct(params.rejectStimuli) && isfield(data,'trial_markings') && size(data.trial_markings,2)>3
 % reject based on stimuli
 items=fieldnames(params.rejectStimuli);
 for i=1:length(items)
  if params.rejectStimuli.(items{i})==1
   % reject this one
   remTrial=find(strcmp(data.trial_markings(:,4),items{i}))';
   if output
    fprintf('Found %d bad trial(s) based on stimulus=%s\n',length(remTrial),items{i})
   end
   badTrials = [badTrials remTrial];
  end
 end
end

if ~isempty(badTrials)
 badTrials =  unique(badTrials);
 goodTrials = goodTrials(~ismember(goodTrials,badTrials));
end

if output
 fprintf('\nTotal: GoodTrials, N=%d / BadTrials, N=%d\n',length(goodTrials),length(badTrials))
end

end

