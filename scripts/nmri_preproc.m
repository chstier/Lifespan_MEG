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
 subject.dws_filt_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_' subject.id '_' subject.exam_id '.mat']);
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
if isfield(subject,'dataset_mapping') && exist(subject.dataset_mapping,'file')
 [ params ] = nmri_get_dataset_params( params, subject.dataset_mapping );
 % load the dataset_id
 load(subject.dataset_mapping,'dataset_id')
else
 dataset_id=[];
end

% check the filtered and downsampled version
if (exist(subject.dws_filt_dataset,'file'))
 load(subject.dws_filt_dataset);   
else
  
 % check what data we have
 switch subject.dtype
  case 'MEG'
   active_channels=ft_channelselection(subject.valid_channels,subject.hdr.grad.label);
  case 'EEG'  
   % Fieldtrip has some trouble with lower/upper cases 
   % make out own channel selection based on header
 
   % try with upper
   active_channels1=ft_channelselection(subject.valid_channels,upper(subject.hdr.label));
   % there may be issues with atypical channel labels, use chantype instead
   active_channels2=subject.hdr.label(strcmpi(subject.hdr.chantype,'EEG'));
   
   % and take the bigger cohort...
   if length(active_channels1)>length(active_channels2)
    active_channels=active_channels1;
   else
    active_channels=active_channels2;
   end
   
   
   % and map back to real
   for i=1:length(active_channels)
    indx=find(strcmpi(subject.hdr.label,active_channels{i}));
    % and make sure that we only have eeg type channels (for special
    % scenarios)
    if ~strcmpi(subject.hdr.chantype{indx},'eeg')
     active_channels{i}=[];
    else
     active_channels{i}=subject.hdr.label{indx};
    end
   end
   % remove empty
   active_channels(cellfun(@isempty,active_channels))=[];
   
   case 'EEG_invasive'
   % Fieldtrip has some trouble with lower/upper cases
   % take upper
   
   % make out own mapping based on chantype
   active_channels=ft_channelselection(strcat(upper(subject.found_electrodes),'*'),upper(subject.hdr.label));
   % and map back to real
   for i=1:length(active_channels)
    indx=find(strcmpi(subject.hdr.label,active_channels{i}));
    % and make sure that we only have eeg type channels (for special
    % scenarios)
    if ~strcmpi(subject.hdr.chantype{indx},'eeg')
     active_channels{i}=[];
    else
     active_channels{i}=subject.hdr.label{indx};
    end
   end
   % remove empty
   active_channels(cellfun(@isempty,active_channels))=[];
   
  otherwise
   error('Unsupported dataset - problem with selecting the active channels')
 end
  
 % check that we have some channels
 if (~exist('active_channels','var') || isempty(active_channels))
  error ('No active channels found in your data / preset combination ')
 end
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

cN=length(active_channels);
cn=1;
chunk_c=0;
if (isfield(params,'mem_max') && ( ~exist('mff_reader','var') || strcmp(mff_reader,'egi_mff_v1'))) % the JAR reader cannot read channel by channel...
 % auto-determine channel_chunks by memory target
 needed_per_channel=subject.hdr.nSamples*subject.hdr.nTrials*32;
 params.channel_chunks=min(floor(params.mem_max/needed_per_channel),subject.hdr.nChans);
 fprintf('Mem need per channel=%.1fM, setting chunk size=%d\n',needed_per_channel/(1024*1024),params.channel_chunks)
elseif (~isfield(params,'channel_chunks')) 
 params.channel_chunks=9999;
end
data_ds=cell(ceil(cN/params.channel_chunks),1); 
%cycle through channels
while (cn<cN)
 chunk_c=chunk_c+1;
 cn_end=cn+params.channel_chunks-1;
 if (cn_end>cN)
  cn_end=cN;
 end
 fprintf('Processing channel(s) %d - %d\n',cn,cn_end)
 %preprocess 
 cfg.channel = active_channels(cn:cn_end);
 data        = ft_preprocessing(cfg);

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

%% now check if we need to remap the channel labels
if isfield(dataset_id,'label_original')
 for i=1:length(data.label)
  if ~any(strcmpi(data.label{i},dataset_id.label)) && sum(strcmpi(data.label{i},dataset_id.label_original))==1
   % not a valid label and one match, remap
   new_i=find(strcmpi(data.label{i},dataset_id.label_original));
   fprintf('Remapping original channel=%s to %s\n',data.label{i},dataset_id.label{i})
   data.label{i}=dataset_id.label{i};
  end
 end
end

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
nmri_write_dataset(subject.dws_filt_dataset,data,subject);

%% deal with headmovement for CTF MEG
if (isfield(params,'track_CTF_movements') && params.track_CTF_movements==1 &&...
  any(strcmpi(data.label,'HLC0011')))
 disp('Now read in coil movements from CTF file')
 cfg=[]; 
 cfg.dataset  = subject.raw_dataset;
 cfg.channel                 = {'HLC0011','HLC0012','HLC0013', 'HLC0021','HLC0022','HLC0023', 'HLC0031','HLC0032','HLC0033'};
 headpos = ft_preprocessing(cfg);

 if (isfield(params,'trial_length') && isnumeric(params.trial_length) && params.trial_length>0) 
  % check if we need to resample
  if abs((size(headpos.trial{1},2)/headpos.fsample)-params.trial_length)>0.01
   cfg           = [];
   cfg.length    = params.trial_length; %single number (in unit of time, typically seconds) of the required snippets
   cfg.overlap   = 0;
   headpos         = ft_redefinetrial(cfg,headpos);
  end
 end
% calculate the mean coil position per trial
 ntrials = size(headpos.sampleinfo,1);
 for t = 1:ntrials
  coil1(:,t) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
  coil2(:,t) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
  coil3(:,t) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
 end

 % calculate the headposition and orientation per trial (for function see bottom page) 
 cc = circumcenter(coil1, coil2, coil3);
 cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';

 % plot translations
 hFig=figure('Position',[0,0,600,450],'Visible','off');
 subplot(2,1,1)
 plot(cc_rel(:,1:3)*1000) % in mm
 title({[subject.id ' ' strrep(subject.exam_id,'_','-')],'Translations in mm'})
 % plot rotations
 subplot(2,1,2) 
 plot(cc_rel(:,4:6))
 title({[subject.id ' ' strrep(subject.exam_id,'_','-')],'Rotations in degree'})
 maxposchange = max(abs(cc_rel(:,1:3)*1000)); % in mm

 subject.movement.cc=cc;
 subject.movement.cc_rel=cc_rel;
 subject.movement.maxposchange=maxposchange;
 if (isfield(params,'QC_plots') && params.QC_plots==1)
  % make the QC dir
  QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
  if (~exist(QCdir,'dir'))
   mkdir(QCdir)
  end
  saveas(hFig,fullfile(QCdir,['movement_CTF_coils_' subject.id '_' subject.exam_id '.png']),'png'); 
 end
 set(hFig,'Visible','on')
end




    

end

