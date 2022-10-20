function [ events, infos ] = nmri_read_event_mff( xmlfile, subject )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~exist(xmlfile,'file')
 error('XMLFile not found')
end

if ~exist('subject','var') || ~isfield(subject,'hdr') || ~isfield(subject,'raw_dataset') || ~exist(subject.raw_dataset,'dir')
 error('No access to raw dataset')
end

% new get the start time from the original dataset
infofile=fullfile(subject.raw_dataset,'info.xml');
if ~exist(infofile,'file')
 error('Could not find info.xml in .mff dataset')
end

xml = nf_xml2struct(infofile);
begTime = xml.recordTime;
begSDV = datetime(begTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX','TimeZone','Europe/Berlin');


% nwo read the xml
xml = nf_xml2struct(xmlfile);


% now read (spike) events (SPK1)
events.evt_timings_seconds={};
events.evt_IDs={};

% and the markings (non spike)
infos.evt_timings_seconds={};
infos.evt_IDs={};

for i=1:length(xml)
 if isfield(xml(i),'event') &&  isfield(xml(i).event,'beginTime')
  eventTime  = xml(i).event.beginTime;
  if ~strcmp('-',eventTime(21)) && isfield(xml(i),'event') && isfield(xml(i).event,'label') && ischar(xml(i).event.label) && ~isempty(xml(i).event.label)
   % event out of range (before recording started): do nothing.
   eventSDV=datetime(eventTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX','TimeZone','Europe/Berlin');
   eventOffset = seconds(eventSDV - begSDV); % duration object in seconds
   % 
   if length(xml(i).event.code)>2 && strcmpi(xml(i).event.code(1:3),'SPK')
    % spike marking
    events.evt_timings_seconds{end+1,1}=eventOffset;
    events.evt_IDs{end+1,1}=xml(i).event.label;
   else
    % something else
    infos.evt_timings_seconds{end+1,1}=eventOffset;
    infos.evt_IDs{end+1,1}=xml(i).event.label;
   end
  end
 end
end
end