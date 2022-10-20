function [subject] = cs_preproc_physio(subject)
% This function does downsampling and trial segmentation for EOG and ECG channels
% Written and modified by Christina Stier, 2020/22

%% Get subject and check whether downsampling of ECG and EOG data has been done already
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

subject2 = subject;
% make our output path and dir
if (~isfield(subject2,'dws_filt_physio'))
 subject2.dws_filt_physio=fullfile(subject2.analysis_dir,subject2.id,'processed',['dws_filt_physio_' subject2.id '_' subject2.exam_id '.mat']);
end
if (~exist(fullfile(subject.analysis_dir,subject.id,'processed'),'dir'))
 mkdir(fullfile(subject.analysis_dir,subject.id,'processed'))
end

% call central dataset selector
subject2=nmri_determine_datatype(subject2);
if isfield(subject2,'dataset_mapping')
 [ params ] = nmri_get_dataset_params( params, subject2.dataset_mapping );
end
 
%% Get the modality-specific analysis params, either via dataset mapping or general and update params
[ params ] = nmri_get_modality_params( params, subject2.dtype );

%% skip makeing the trials here, can do later
cfg          = params.preproc_cfg;
cfg.dataset  = subject2.raw_dataset;

if strcmpi(subject2.raw_dataset(end-3:end),'.mff') || ( isfield(subject2,'detected_datatype') && strcmp(subject2.detected_datatype(1:4),'EGI-'))
 % use v2/v3 mff reader for EGI here, clipped files are not read by v1...so
 % try to be smart with the central function
 mff_reader=nmri_check_mff_reader(subject2.raw_dataset);
 fprintf('Using MFF-Reader=%s\n',mff_reader);
 cfg.headerformat=mff_reader;
 cfg.dataformat=mff_reader;
 subject2.mff_reader=mff_reader;
end

% unless we have a specific trial function ... then do now
if (isfield(params,'trial_fun') && ischar(params.trial_fun) && ~isempty(params.trial_fun))
 cfg.trialfun = params.trial_fun; % ft_definetrial will call your function and pass on the cfg
 cfg          = ft_definetrial(cfg);
end

if (isfield(params,'dws_freq')) 
 ds_cfg            = [];
 ds_cfg.detrend    = 'no'; % generally not used
 ds_cfg.resamplefs = params.dws_freq;
end

if (isfield(params,'dws_freq')) 
 ds_cfg            = [];
 ds_cfg.detrend    = 'no'; % generally not used
 ds_cfg.resamplefs = params.dws_freq;
end

 cfg.channel = {'eog','ecg'};
 data2        = ft_preprocessing(cfg);
 
 data_ds = ft_resampledata(ds_cfg, data2);
 data2 = data_ds;
  % make sure we keep the hdr
 if isfield(data_ds,'hdr')
  data2.hdr=data_ds.hdr;
 end
 
 
 %% downsampling destroys sampleinfo - regenerate
w=warning('off','all');
data2=ft_checkdata(data2, 'feedback', 'no', 'hassampleinfo', 'yes');
warning(w);

%% now we redefine our trials

if (isfield(params,'trial_length') && isnumeric(params.trial_length) && params.trial_length>0) 
 cfg           = [];
 cfg.length    = params.trial_length; %single number (in unit of time, typically seconds) of the required snippets
 cfg.overlap   = 0;
 data2         = ft_redefinetrial(cfg,data2);
end

%% make all channels upper case - we not need to refer to the original datset any more
data2.label=upper(data2.label);
subject2.hdr.label=upper(subject2.hdr.label);

%% if we got this far, place a stamp for completed pre-processing
subject2=nmri_stamp_subject(subject2,'preproc_physio',params);
subject=nmri_stamp_subject(subject,'preproc_physio',params);

%% save downsampled data and subject
nmri_write_dataset(subject2.dws_filt_physio,data2,subject2);

%% visual check
ft_databrowser([],data2);