function [ subject, ekg ] = nmri_read_ekg(subject, params)
%[ subject, ekg ] = nmri_read_ekg(subject, params)
%  
% This function will try read an EKG/ECG channel from MEG/EEG datasets
% in a way that is compatiple to nmri_preproc
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% params    =   anaylsis parameter struct (optional, will search for
%               analysis_params.m if not set)

% written by NF 09/2019

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
 subject.dws_filt_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_' subject.id '_' subject.exam_id '.mat']);
end
if (~isfield(subject,'ekg_dataset'))
 subject.ekg_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['EKG_' subject.id '_' subject.exam_id '.mat']);
end
if (~exist(fullfile(subject.analysis_dir,subject.id,'processed'),'dir'))
 mkdir(fullfile(subject.analysis_dir,subject.id,'processed'))
end

if (~exist(subject.dws_filt_dataset,'file'))
 error('Could not find dws_filt dataset. Needed to time-sync ECG');   
else


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
 % use v2 mff reader for EGI here, clipped files are not read by v1...so
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
 
%% now do basic preprocessing

try
 cfg.channel = {'ECG','EKG'}; 
 fprintf('\nReading EKG...\n')
 ekg        = ft_preprocessing(cfg);
catch ME
 if strcmpi(ME.message,'the selection of channels is empty')
  error('No ECG/EKG channel was found in the dataset by Fieldtrip, check that this channel really exists')
 else
  rethrow(ME)
 end
end

% note, there is a problem in Fieldtrip that with the mff_v3 reader all
% PNS channels are concatenated with the EEG without correct referencing,
% FIX: if so, remap to the original channels from egi_mff_v3
if isfield(subject,'mff_reader') && strcmpi(subject.mff_reader,'egi_mff_v3')...
  && size(ekg.trial{1},1)~=size(ekg.label,1)
 fprintf('EGIv3 reader has mis-matched the channel labels, remapping from header\n')
 ekg.label=subject.hdr.label;
 % re-check
 if size(ekg.trial{1},1)~=size(ekg.label,1)
  error('Could not fix the divergent number of channels read from the EGI file, investigate here...')
 end
end

% now make sure to select EKG only
if size(ekg.label,1)>1
 cfg=[];
 cfg.channel = {'ECG','EKG'}; 
 ekg=ft_selectdata(cfg,ekg);
end

% resample
if (isfield(params,'dws_freq')) 
 ds_cfg            = [];
 ds_cfg.detrend    = 'no'; % generally not used
 ds_cfg.resamplefs = params.dws_freq;
 ekg = ft_resampledata(ds_cfg, ekg);
end

%% downsampling destroys sampleinfo - regenerate
w=warning('off','all');
ekg=ft_checkdata(ekg, 'feedback', 'no', 'hassampleinfo', 'yes');
warning(w);

%% now we redefine our trials

if (isfield(params,'trial_length') && isnumeric(params.trial_length) && params.trial_length>0) 
 cfg           = [];
 cfg.length    = params.trial_length; %single number (in unit of time, typically seconds) of the required snippets
 cfg.overlap   = 0;
 ekg         = ft_redefinetrial(cfg,ekg);
end

%% make all channels upper case - we not need to refer to the original datset any more
ekg.label=upper(ekg.label);

% and make sure we have only one
if length(ekg.label)~=1
 cfg=[];
 % take the first
 cfg.channel=ekg.label{1};
  fprintf('\nFound more then one match, w EKG...\n')
 ekg=ft_selectdata(cfg,ekg);
end


% and make sure we have a compatible time to ds_filt
load(subject.dws_filt_dataset,'data'); 
if ~isequal(ekg.time,ekg.time)
 disp('Timing not equal beteen EKG and data, resampling EKG to match data')
 cfg=[];
 cfg.time=ekg.time;
 ekg=ft_resampledata(cfg, ekg);
end

%% if we got this far, place a stamp for completed EKG reading
subject=nmri_stamp_subject(subject,'ekg-read',params);

%% save downsampled data and subject

save(subject.ekg_dataset,'ekg','subject');



    

end

