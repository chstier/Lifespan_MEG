function nmri_write_dataset(fname, data, subject)
%nmri_write_dataset(fname, data, subject)
%  Function to write a dataset to disk as .mat file
% we also generate some "smart" separate structs for fast access
% and deal with synchonising subject and data struct
%
% fname = filename
% data = fieldtrip data structure
% subject = subject structure
% ... further variables may also be provided and are saved as is
%

if ~exist('fname','var') || ~exist('data','var') || ~exist('subject','var')
 error('Always need filename, data and subject')
end

data_info=[];
% check if trial markings are present and save
if isfield(data,'trial_markings')
 data_info.trial_markings=data.trial_markings;
end
if isfield(data,'trial_markings_sampleinfo')
 data_info.trial_markings_sampleinfo=data.trial_markings_sampleinfo;
end
if isfield(data,'bad_channels')
 data_info.bad_channels=data.bad_channels;
end

data_info.ntrials=size(data.trial,2);
data_info.nchannels=size(data.trial{1},1);
data_info.fsample=data.fsample;

subject.data_info=data_info;
fprintf('Saving to file=%s\n',fname)
save(fname,'data','subject')

end

