function  [ events, infos ] = nmri_read_event_edf( edffile, subject )
% [ events, infos ] = nmri_read_event_edf( edffile, subject )
%   Reads in EDF annotations

if ~exist(edffile,'file')
 error('EDF file not found')
end

if ~exist('subject','var') || ~isfield(subject,'hdr') || ~isfield(subject.hdr,'orig') || ~isfield(subject.hdr.orig,'annotation')
 error ('Subject not set or HDR not found or not EDF(+) with annotations')
end

all_evts=ft_read_event(subject.raw_dataset,'header',subject.hdr,'detectflank',[]);

% get the offset (Fieldtrip does nto seem to comply with EDF+ rules in
% 20180801... the first event is the offset

if length(all_evts)>0
 offset=all_evts(1).timestamp;
else
 offset=0;
end

% now read (spike) events (SPK1)
events.evt_timings_seconds={};
events.evt_IDs={};

% and the markings (non spike)
infos.evt_timings_seconds={};
infos.evt_IDs={};

% now parse the events
for i=1:length(all_evts)
 if ~isempty(all_evts(i).value)
  if ~isempty(regexpi(all_evts(i).value,'sp[ike]+'))
   % probably a spike marking
   events.evt_timings_seconds{end+1,1}=all_evts(i).timestamp-offset;
   events.evt_IDs{end+1,1}=all_evts(i).value;
  else
   % something else
   infos.evt_timings_seconds{end+1,1}=all_evts(i).timestamp-offset;
   infos.evt_IDs{end+1,1}=all_evts(i).value;
  end
 end
end

end

