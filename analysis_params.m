%analysis parameter definitions, to be used for all subjects in this stream

% This is a unified configuaration file for EEG (EGI 256) und MEG 
params=[];

% Directory that has the mapping functions for dataset detection
params.dataset_mapping_dir='conf/dataset-mapping'; 


% ----------------------------------------------------------------------
% Preprocessing Options

% Preoprocessing option - as in Fieldtrip ft_preprocessing
params.preproc_cfg.continuous = 'yes'; % we have continuous data for resting-state

% Frequency Filter 
params.preproc_cfg.hpfilter   = 'yes';  % highpass
params.preproc_cfg.hpfreq     = 1;      % Hz
params.preproc_cfg.hpfiltord  = 1; % first order Butterworth
%params.preproc_cfg.lpfilter   = 'yes';  % lowpass 
%params.preproc_cfg.lpfreq     = 70;     % Hz
%params.preproc_cfg.lpfiltord  = 1; % first order Butterworth

% Line noise
params.preproc_cfg.dftfilter  = 'no';   % use (sharp) line noise removal
params.preproc_cfg.dftfreq    = [50 100 150]; % default to use 3 harmonics

% Bandstop filters - seems needed for EGI EEG
params.EEG.preproc_cfg.bsfilter  = 'yes';   % additional line noise removal
params.EEG.preproc_cfg.bsfreq    = [45 55]; % 

% Trend / Mean removal
params.preproc_cfg.detrend    = 'no';
params.preproc_cfg.demean     = 'yes'; % substract the mean, vital for EEG,
% probably no harm done in MEG

% Trial definiation
params.trial_length = 10; % in seconds
%params.trial_fun='sub_trialfun'; % use a specific trial function

% Downsampling
params.dws_freq = 300; % will be sampled to this frequency
% if not set, will keep original sampling rate

% special freq bands for invasive EEG
params.EEG_invasive.preproc_cfg.lpfreq  = 1000;
params.EEG_invasive.dws_freq = 2000; % will be sampled to this frequency

% how many channels to be done in parallel
% depeneds on how much memory you have, if not set, all channels will be
% done at once (as if set to Inf)
%params.channel_chunks = 15; % 1 - 99999 

% or set via auto
params.mem_max=4*1024*1024*1024; % 4 GB


% ----------------------------------------------------------------------
% Frequency and nTrials

% define frequency bands
params.freqs = [2 6 10 16 25 40 ]; % in Hz, center frequency of interest
params.tapsmofrq = [2 2 2 4 4 8 ]; % plus minus the frequency smoothing

params.freqsNames = {'Delta' 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Gamma'};

% select how many trials to use, need to have at least this many in all
% datasets (after rejection), otherwise will fail
params.nTrials =30;



% ----------------------------------------------------------------------
% Artifact cleaning options

% State if you want automated artifact rejection
params.useAUTO_clean = 1; % commment/unset if not wanted or set to 0

% Options for auto clean

% 1st pass - Trials
% params.AUTO_clean.params=[];
% params.AUTO_clean.params(1).mode = 'trials'; % type of rejection
% params.AUTO_clean.params(1).metrics = {'var','var','var','kurtosis','abszvalue'}; % metric and order for rejection
% params.AUTO_clean.params(1).cutoff = [3 2 2 3 3]; % MAD (mean absolute deviation) to use as cutoff, one for each metric

% % 2nd pass - Trials by variance in channels
% params.AUTO_clean.params(2).mode = 'trials_by_channels'; % type of rejection
% params.AUTO_clean.params(2).metrics = {'var','var','var','kurtosis'}; % metric and order for rejection
% params.AUTO_clean.params(2).cutoff = [3 2 1.5 2]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
% 
% % 3rd pass - Channels
% params.AUTO_clean.params(3).mode = 'channels'; % type of rejection
% params.AUTO_clean.params(3).metrics = {'var','var','kurtosis','1/var','1/kurtosis','abszvalue'}; % metric and order for rejection
% params.AUTO_clean.params(3).cutoff = [4 3 4 4 4 4]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
% 
% % 4th pass - Trials again
% params.AUTO_clean.params(4).mode = 'trials'; % type of rejection
% params.AUTO_clean.params(4).metrics = {'kurtosis','abszvalue','var'}; % metric and order for rejection
% params.AUTO_clean.params(4).cutoff = [2 2 2]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
% 
% % 5th pass - Channels again
% params.AUTO_clean.params(5).mode = 'channels'; % type of rejection
% params.AUTO_clean.params(5).metrics = {'var'}; % metric and order for rejection
% params.AUTO_clean.params(5).cutoff = [3]; % MAD (mean absolute deviation) to use as cutoff, one for each metric



params.AUTO_clean.params=[];
params.AUTO_clean.params(1).mode = 'trials'; % type of rejection
params.AUTO_clean.params(1).metrics = {'var','var','var'}; % metric and order for rejection
params.AUTO_clean.params(1).cutoff = [3 3 3]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
params.AUTO_clean.params(2).mode = 'channels'; % type of rejection
params.AUTO_clean.params(2).metrics = {'var','var','var'}; % metric and order for rejection
params.AUTO_clean.params(2).cutoff = [4 3 3]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
params.AUTO_clean.params(3).mode = 'trials'; % type of rejection
params.AUTO_clean.params(3).metrics = {'var','var','var'}; % metric and order for rejection
params.AUTO_clean.params(3).cutoff = [3 2 2]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
params.AUTO_clean.params(4).mode = 'channels'; % type of rejection
params.AUTO_clean.params(4).metrics = {'var','var','var'}; % metric and order for rejection
params.AUTO_clean.params(4).cutoff = [4 3 3]; % MAD (mean absolute deviation) to use as cutoff, one for each metric
params.AUTO_clean.params(5).mode = 'trials'; % type of rejection
params.AUTO_clean.params(5).metrics = {'var','var','var'}; % metric and order for rejection
params.AUTO_clean.params(5).cutoff = [3 2 2]; % MAD (mean absolute deviation) to use as cutoff, one for each metric



% State if you want visual artifact cleaning 
params.useVIZ_clean = 1; % commment/unset if not wanted or set to 0

% State if you want ICA cleaning as well
params.useICA_clean = 1; % commment/unset if not wanted or set to 0



% ----------------------------------------------------------------------
% Events and trial rejection

% Set if you want to have trials rejected that contain events (e.g. spikes
% or special artifacts)
params.rejectEvents = 1;

% Rejection will be based on the dataset mapping in conf/dataset_mappings
% this will override the functions given here!!
params.events_per_dataset = 1; 

% path to marking file
params.MEG.rejectEventsPath = '<analysis_dir>/spike_markings';
% pattern matching, <subject_id> will be replaced by actual subject ID
params.MEG.rejectEventsPattern = '<subject_id>*.evt';
% textscan filter to read the spike time from file
params.MEG.rejectEventsTextscan='[read_cell]=textscan(fid,''%d %s %*[^\n]'',''HeaderLines'',1);spikeTimes=read_cell{1};spikeIDs=read_cell{2}'; % BESA for MEG
% this is a Matlab eval, so can be customized as needed. Has to provide
% spikeTimes, spikeIDs is optional (usefull if more than one spike type)
% timebase of markings, e.g. 1000000 (microseconds for BESA, Tmu)
params.MEG.rejectEventsTimebase = 1000000;


% different for EEG
%EGI-MFF
% path to marking file 
% params.EEG.rejectEventsPath = {'<raw_dataset>','<analysis_dir>/spike_markings'};
% params.EEG.rejectEventsPattern = {'Event*.xml','<subject_id>*.xml'};
% Call a function to read events, 1st parameter is the file, second subject
% params.EEG.rejectEventsFunction='nmri_read_event_mff'; % EGI for EEG
% For EGI EEG already seconds
% params.EEG.rejectEventsTimebase = 1; %


% Use Geoscan or not
params.EEG.Geoscan_Use = 1; % will use the Geoscan, if available
params.EEG.Geoscan_Require = 0; % will require the Geoscan and fail otherwise



%EDF
% path to marking file 
params.EEG.rejectEventsPath = {'<raw_dataset>'};
params.EEG.rejectEventsPattern = {'.'};
% Call a function to read events, 1st parameter is the file, second subject
params.EEG.rejectEventsFunction='nmri_read_event_edf'; % 
% For EEG already seconds
params.EEG.rejectEventsTimebase = 1; %


% Include an automated mapping of stimuli during read_markers
params.parse_stimuli = 1; % 1: do it, 0: skip it 


% trials to exclude, relative to affected
params.rejectEventsTrial = [-1 0 1]; % will remove trial before and after
% will map based on trail length/sampleinfo

% period for spike analysis (PrePost)
params.EventsTOI = 5; % in s

% how to deal with vigilance scoring and trial selection
params.scoreVigilance = 1; % request to score vigilance in all subjects

% which vigilance levels to accept
params.rejectVigilance.wake = 0; % only accept wake
params.rejectVigilance.sleep1 = 1;
params.rejectVigilance.sleep2 = 1;
params.rejectVigilance.sleep3 = 1;
params.rejectVigilance.rem = 1;


% which stimulation levels to accept
% params.rejectStimuli.eyes_closed = 0; % accept eyes-closed
% params.rejectStimuli.none = 0; % accept no score
% params.rejectStimuli.eyes_open = 1; % remove eyes-open
% params.rejectStimuli.HV = 1; % remove HV
% params.rejectStimuli.postHV = 1; % remove postHV
% params.rejectStimuli.PS = 1; % remove photo-stim


% Default electrode layout scheme (subject may overwrite)
% now detected by nmri_determine_datatype and mapping function (above)
%params.MEG.layout   = 'CTF275.lay'; %MEG default
%params.EEG.layout   = 'conf/GSN-HydroCel-257-layout.mat'; % EGI 256 (+Ref) default
%params.EEG.layout   = 'conf/EEG1020_Goe.lay'; % 10-20 default, for routine EEG

% Electrodes file, only needed for EEG
%params.EEG.elec_file   = 'GSN-HydroCel-257.sfp';
% now also give fiducial points of this config - Corresponding to NAS, LPA
% and RPA (in this order)
%params.EEG.elec_fiducials = {'FidNz','FidT9','FidT10'}; % for EGI 256 (+Ref) cap



% ----------------------------------------------------------------------
% SUMA / Headmodel

% SUMA surface to use
params.SUMA_ld='10';
params.SUMA_surf='smoothwm';

% Segmentation controls - use multichannel (FLAIR or T2)?
% if so, these need to be run through Freesurfer as well (deals with bias
% correction and registration already)
params.use_FLAIR = 0; % will use FLAIR as well
params.use_FLAIR_permissive = 1; % but will not fail, if no FLAIR is found

params.use_T2 = 1; % will use T2 as well - change for CamCan datasets
params.use_T2_permissive = 1; % but will not fail, if no T2 is found 

% Also use subcortical nuclei? 
params.ASEG_use=1;
% Do you want fully individual subcortical ROIs?
% otherwise will do MNI-space mapping of fsaverage ROIs
% MNI-mapped has the advantage to be vertex-wise comparable between
% subjects but could have anantomical imprecissions.  
params.ASEG_individual = 0; % 1= individual, 0=MNI-mapped ROIs
% requires a configuation text file with the SUMA ld as above
params.ASEG_conf=['/tools/common/freesurfer-labels-aseg-connectome-' num2str(params.SUMA_ld)];
% shall we use the "old style" ASEG / subcortial matching based on SPM
% unified segmenation? (default as of 20170418 is to use DARTEL instead)
% for the DARTEL stream to work, we need the DARTEL templates in
% /tools/common/DARTEL_templates_CAT12
params.ASEG_MNI_unified=0; % set to 1 for old style/unified segmentation

% MRI parameters
params.MRI.dim=[256 256 256]; % dimensions (voxels of MRI)
% default for Freesurfer is [256 256 256] for 1mm voxels
% but this seems to small in some cases
% if you use high-res data this needs to be increased further

% headmodel
% for MEG
params.MEG.headmodel = 'singleshell'; % usual MEG - "Nolte" method
params.MEG.headmodel_class = 'brain'; % 'brain' or 'suma'
% usual behaviour of FT is to use 'brain', but 'suma' seems to be very
% similar (in this way the SUMA surface is used directly as solution
% points, i.e. the point are ON the singleshell, not INSIDE)

% for EEG
params.EEG.headmodel = 'dipoli'; % BEM headmodel
%params.EEG.headmodel = 'simbio'; % FEM headmodel

% cutoff for tissue class generation from SPM c5 class
params.c5cutoff = 0.9; %0.5-0.99 would be reasonable
% lower values can cover "holes" in the MRI with the information in the
% priors. But can also cause unreal tissue in the headmodel.
% Thus, you may need to modify this parameter depending on the scenario
% Should be safe to modify on a subject-by-subject basis

% for EEG
params.EEG_invasive.headmodel = 'dipoli'; % BEM headmodel
%params.EEG_invasive.headmodel = 'simbio'; % FEM headmod


% ----------------------------------------------------------------------
% Connectivity options

% Connectivity measures to calculate
params.con_method={'wpli_debiased','coh'};

% safe compact stats file (without filters and Fourier
params.stats_compact = 0; % default, safe full file
%params.stats_compact = 1; % safe compact file



% ----------------------------------------------------------------------
% QC options

% make QC plots
params.QC_plots = 1; 

% Extract headmovement parameters - for MEG
params.MEG.track_CTF_movements = 1; % will read CTF coil tracking



% ----------------------------------------------------------------------
% Source projection options

% Options for source projection of timecourses
params.src_project.beamformer='DICS';
%params.src_project.beamformer='LCMV'; %optional







