function [ subject, data ] = nmri_preproc(subject, params)
%[ subject, data ] = nmri_preproc(subject, params)
%  
% This function will do basic preprocessing on one subject
% i.e. trial cutting, filtering and downsampling 
% Usually this is the first step of a processing pipeline
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% params    =   anaylsis parameter struct (optional, will search for
%               analysis_params.m if not set)
% data      =   will return the data (optional)

% written by NF 11/2016 - 09/2018

% this functions need the signal processing toolbox for resampling


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

% check if we have raw_dataset (this is the minimum)
if (~isfield(subject,'raw_dataset') || ~exist(subject.raw_dataset,'file'))
 error('Raw dataset not specified or not accessible, check nmri_read_subject')
end

% make our output path and dir
if (~isfield(subject,'dws_filt_dataset'))
 subject.dws_filt_dataset_refchannels=fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_' subject.id '_' subject.exam_id '.mat']);
end
if (~exist(fullfile(subject.analysis_dir,subject.id,'processed'),'dir'))
 mkdir(fullfile(subject.analysis_dir,subject.id,'processed'))
end

% % read hdr, if not there
% if (~isfield(subject,'hdr') || isempty(subject.hdr))
%  if strcmpi(subject.raw_dataset(end-3:end),'.mff')
%   % use v2 mff reader for EGI here, clipped files are not read by v1...
%   % but old files nor not read by the JAR version, so determine 
%   subject.hdr=ft_read_header(subject.raw_dataset,'headerformat',nmri_check_mff_reader(subject.raw_dataset));
%   % need to reomve java Objects...these do not serialize
%   if isfield(subject,'hdr') && isfield(subject.hdr,'orig') && isfield(subject.hdr.orig,'javaObjs')   
%    subject.hdr.orig=rmfield(subject.hdr.orig,'javaObjs');
%   end
%   % also get rid of MFF_v3 original data...stupid idea to store that here
%   if isfield(subject,'hdr') && isfield(subject.hdr,'orig') && isfield(subject.hdr.orig,'data')   
%    subject.hdr.orig=rmfield(subject.hdr.orig,'data');
%   end
%  else
%   subject.hdr=ft_read_header(subject.raw_dataset);
%  end
% end

% call central dataset selector
subject=nmri_determine_datatype(subject);
if isfield(subject,'dataset_mapping')
 [ params ] = nmri_get_dataset_params( params, subject.dataset_mapping );
end
 
 

 %% Get the modality-specific analysis params, either via dataset mapping or general
 % and update params
[ params ] = nmri_get_modality_params( params, subject.dtype );

%% skip makeing the trials here, can do later
cfg          = params.preproc_cfg;
cfg.dataset  = subject.raw_dataset;

if strcmpi(subject.raw_dataset(end-3:end),'.mff') || ( isfield(subject,'detected_datatype') && strcmp(subject.detected_datatype(1:4),'EGI-'))
 % use v2/v3 mff reader for EGI here, clipped files are not read by v1...so
 % try to be smart with the central function
 mff_reader=nmri_check_mff_reader(subject.raw_dataset);
 fprintf('Using MFF-Reader=%s\n',mff_reader);
 cfg.headerformat=mff_reader;
 cfg.dataformat=mff_reader;
 subject.mff_reader=mff_reader;
end

% unless we have a specific trial function ... then do now
if (isfield(params,'trial_fun') && ischar(params.trial_fun) && ~isempty(params.trial_fun))
 cfg.trialfun = params.trial_fun; % ft_definetrial will call your function and pass on the cfg
 cfg          = ft_definetrial(cfg);
end
 
%% now do basic preprocessing per channel (to save memory)
if (isfield(params,'dws_freq')) 
 ds_cfg            = [];
 ds_cfg.detrend    = 'no'; % generally not used
 ds_cfg.resamplefs = params.dws_freq;
end

% cN=length(active_channels);
% cn=1;
% chunk_c=0;
% if (isfield(params,'mem_max') && ( ~exist('mff_reader','var') || strcmp(mff_reader,'egi_mff_v1'))) % the JAR reader cannot read channel by channel...
%  % auto-determine channel_chunks by memory target
%  needed_per_channel=subject.hdr.nSamples*subject.hdr.nTrials*32;
%  params.channel_chunks=min(floor(params.mem_max/needed_per_channel),subject.hdr.nChans);
%  fprintf('Mem need per channel=%.1fM, setting chunk size=%d\n',needed_per_channel/(1024*1024),params.channel_chunks)
% elseif (~isfield(params,'channel_chunks')) 
%  params.channel_chunks=9999;
% end
% data_ds=cell(ceil(cN/params.channel_chunks),1); 
% %cycle through channels
% while (cn<cN)
%  chunk_c=chunk_c+1;
%  cn_end=cn+params.channel_chunks-1;
%  if (cn_end>cN)
%   cn_end=cN;
%  end
%  fprintf('Processing channel(s) %d - %d\n',cn,cn_end)
 %preprocess 
%  cfg.channel = active_channels(cn:cn_end);
 cfg.channel = {'eog','ecg'};
 data2        = ft_preprocessing(cfg);

 % note, there is a problem in Fieldtrip that with the mff_v3 reader all
 % PNS channels are concatenated with the EEG without correct referencing,
 % FIX: if so, remap to the original channels from egi_mff_v3
 if isfield(subject,'mff_reader') && strcmpi(subject.mff_reader,'egi_mff_v3')...
   && size(data.trial{1},1)~=size(data.label,1)
  fprintf('EGIv3 reader has mis-matched the channel labels, remapping from header\n')
  data.label=subject.hdr.label;
  % re-check
  if size(data.trial{1},1)~=size(data.label,1)
   error('Could not fix the divergent number of channels read from the EGI file, investigate here...')
  end
 end
 
 % now make sure to select valid channels only
 if size(data.label,1)~=length(active_channels)
  scfg=[];
  scfg.channel = active_channels; 
  data=ft_selectdata(scfg,data);
 end

 
 
 if (isfield(params,'dws_freq')) 
  %now downsample
  data_ds{chunk_c} = ft_resampledata(ds_cfg, data);
 else
  data_ds{chunk_c} = data;
 end

 cn=cn+params.channel_chunks;
end
if (chunk_c==1)
 data=data_ds{1};
else
 % more than one chunck, append channels then
 data=ft_appenddata([],data_ds{:});
 % make sure we keep the hdr
 if isfield(data_ds{1},'hdr')
  data.hdr=data_ds{1}.hdr;
 end
end

%% downsampling destroys sampleinfo - regenerate
w=warning('off','all');
data=ft_checkdata(data, 'feedback', 'no', 'hassampleinfo', 'yes');
warning(w);

%% now we redefine our trials

if (isfield(params,'trial_length') && isnumeric(params.trial_length) && params.trial_length>0) 
 cfg           = [];
 cfg.length    = params.trial_length; %single number (in unit of time, typically seconds) of the required snippets
 cfg.overlap   = 0;
 data         = ft_redefinetrial(cfg,data);
end

%% make all channels upper case - we not need to refer to the original datset any more
data.label=upper(data.label);
subject.hdr.label=upper(subject.hdr.label);

%% if we got this far, place a stamp for completed pre-processing
subject=nmri_stamp_subject(subject,'preproc',params);

%% save downsampled data and subject
nmri_write_dataset(subject.dws_filt_dataset_refchannels,data,subject);





    

end

