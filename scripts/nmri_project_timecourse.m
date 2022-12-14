function [ subject ] = nmri_project_timecourse(subject, data, params)
%[ subject ] = nmri_project_timecourse(subject, data, params)
%  
% This function will do the main processing including
% power @sensor level (for QC)
% source recon
% power @source
% connectivity @source
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 11/2016 - 03/2017


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
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


% we do not need trial markings any more now, remove to avoid Fieldtrip
% warnings
if isfield(data,'trial_markings')
 trial_markings=data.trial_markings;
 data=rmfield(data,'trial_markings');
end
if isfield(data,'trial_markings_sampleinfo')
 trial_markings_sampleinfo=data.trial_markings_sampleinfo;
 data=rmfield(data,'trial_markings_sampleinfo');
end

if (strcmp(subject.dtype,'EEG'))
 if (~isfield(subject,'electrodes_aligned') || ~exist(subject.electrodes_aligned,'file'))
  error('Have not found EEG electrodes - should be generated by the headmodel step')
 end
 load(subject.electrodes_aligned,'elec_aligned','elec_present','elec_missing')
end

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% make dir / files
if (~isfield(subject,'source_dir'))
 subject.source_dir=fullfile(subject.analysis_dir,subject.id,'processed','source_projected');
end
if (~exist(subject.source_dir,'dir'))
 mkdir(subject.source_dir)
end

%% load the head model and leadfield
if ~isfield(subject,'hdm_lead') || ~exist(subject.hdm_lead,'file')
 error('Could not find headmodel - Make sure to run nmri_make_hdm_suma first')
end
if ~isfield(subject,'suma_surface') || ~exist(subject.suma_surface,'file')
 error('Could not find SUMA surface - Make sure to run nmri_make_hdm_suma first')
end

load(subject.hdm_lead,'hdm','leadfield');
load(subject.suma_surface,'suma_all');

% check SUMA and leadfield N's
if size(suma_all.pos,1)~=sum(leadfield.inside)
 warning('Not all SUMA points are inside the brain, there may be issues with vertex-dipole matching - investigate')
end

%% Now do the source projection
switch params.src_project.beamformer
 case 'DICS'
  % now loop per frequency for DICS
  for ff = 1:length(params.freqs)
   % create results file
   result_file = fullfile(subject.source_dir,[params.src_project.beamformer '_' params.freqsNames{ff} '.mat']);
   if ~exist(result_file,'file')
    % not present, then calculate  
    fprintf('DICS: Working on Frequency=%s\n',params.freqsNames{ff})

    cfg           = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.pad       = 'nextpow2'; % for more efficient FFT computation
    cfg.foi       = params.freqs(ff);
    cfg.tapsmofrq = params.tapsmofrq(ff); 
    tmpfreq       = ft_freqanalysis(cfg, data);
    
    cfg                 = [];
    cfg.channel         = data.label;
    cfg.method          = 'dics';
    cfg.pad             = 'nextpow2'; % recommended for speed
    cfg.frequency       = params.freqs(ff);
    cfg.headmodel       = hdm;
    cfg.grid            = leadfield;
    cfg.dics.lambda     = '5%';
    cfg.dics.keepfilter = 'yes';
    cfg.dics.fixedori   = 'yes';
    cfg.dics.realfilter = 'yes';
    cfg.dics.projectnoise  = 'yes';
    cfg.dics.normalize  = 'yes';
    
    % for EEG need the electrodes
    if (strcmp(subject.dtype,'EEG'))
     cfg.elec            = elec_present;
    end
    
    tmpsource           = ft_sourceanalysis(cfg, tmpfreq);
    
    data_source=data;
    
    % make labels for SUMA sources
    data_source.label = cellstr(num2str((1:sum(leadfield.inside))'));
        
    % now project the timecourse into source space - NOTE: these are
    % only INSIDE dipoles. i.e there may be a mismatch to SUMA in rare cases
   
    fprintf('DICS: Projecting data now\n')
    
    fprintf('Source: ')
    ss=0;
    for i = 1:length(tmpsource.avg.filter) % loop for each source/filter
     if i>1
      for j=0:log10(i-1)
       fprintf('\b',i)
      end
     end
     fprintf('%d',i)
     if ~isempty(tmpsource.avg.filter{i})
      % for non-inside the filter will be empty, but still present
      ss=ss+1;
      % now loop for trials
      for iTrial = 1:length(data.trial)
       data_source.trial{iTrial}(ss,:) = tmpsource.avg.filter{ss} * data.trial{iTrial} ;
      end
     end
    end
    fprintf('...done\n\n')
    
    if (ss~=sum(leadfield.inside))
     error('Have not found the needed N of inside sources. This should not happen')
    end
    
    % map to SUMA if needed
    if size(suma_all.pos,1)~=ss
     warning('Mismatch of source-points and SUMA, will try to re-map')
     for iTrial = 1:length(data_source.trial)      
      mapped_trial=zeros(size(suma_all.pos,1),size(data_source.trial{iTrial},2));
      mapped_trial(leadfield.inside,:)=data_source.trial{iTrial};
      mapped_trial(~leadfield.inside,:)=NaN;
      mapped_trial(:,~leadfield.inside)=NaN;
      data_source.trial{iTrial}=mapped_trial;
     end
    end
    
    % get back trial markings
    if exist('trial_markings_sampleinfo','var')
     data_source.trial_markings_sampleinfo=trial_markings_sampleinfo;
    end
    if exist('trial_markings','var')
     data_source.trial_markings=trial_markings;
    end
    
    % save data
    save(result_file,'data_source','subject')
   else
    fprintf('DICS: Data present for Frequency=%s\n',params.freqsNames{ff})
   end
  end
  
 case 'LCMV'
  % create results file
  result_file = fullfile(subject.source_dir,[params.src_project.beamformer '.mat']);
  if ~exist(result_file,'file')
   % not present, then calculate  
   fprintf('LCMV: Generating broadband\n')
   cfg              = [];
   cfg.covariance   ='yes';
   cfg.vartrllength = 2;% accept variable trial length
   avg = ft_timelockanalysis(cfg,data);
 
   cfg                     = [];
   cfg.method              = 'lcmv';
   cfg.grid                = leadfield;
   cfg.headmodel           = hdm;
   cfg.lcmv.keepfilter     = 'yes';
   cfg.lcmv.lambda         = '5%';
   cfg.lcmv.weightnorm     = 'nai';  % this is a new option :-) from Saran Dalal in mailing list
   cfg.lcmv.fixedori       = 'yes';
   cfg.lcmv.projectnoise   = 'yes';
   cfg.lcmv.projectmom     = 'no';
   cfg.lcmv.reducerank     = 'no';
   cfg.channel             = data.label;
   
   % for EEG need the electrodes
   if (strcmp(subject.dtype,'EEG'))
    cfg.elec            = elec_present;
   end
   
   tmpsource               = ft_sourceanalysis(cfg, avg);
  
   data_source=data;
   
   % make labels for SUMA sources
   data_source.label = cellstr(num2str((1:sum(leadfield.inside))'));
        
   % now project the trial timecourse into source space - NOTE: these are
   % only INSIDE dipoles. i.e there may be a mismatch to SUMA in rare cases
   fprintf('LCMV: Projecting data now\n')
   fprintf('Source: ')
   ss=0;
   for i = 1:length(tmpsource.avg.filter) % loop for each source/filter
    if i>1
     for j=0:log10(i-1)
      fprintf('\b',i)
     end
    end
    fprintf('%d',i)
    if ~isempty(tmpsource.avg.filter{i})
     % for non-inside the filter will be empty, but still present
     ss=ss+1;
     % now loop for trials
     for iTrial = 1:length(data.trial)
      data_source.trial{iTrial}(ss,:) = tmpsource.avg.filter{ss} * data.trial{iTrial} ;
     end
    end
   end
   fprintf('...done\n\n')

   if (ss~=sum(leadfield.inside))
    error('Have not found the needed N of inside sources. This should not happen')
   end
    
   % map to SUMA if needed
   if size(suma_all.pos,1)~=ss
    warning('Mismatch of source-points and SUMA, will try to re-map')
    for iTrial = 1:length(data_source.trial)      
     mapped_trial=zeros(size(suma_all.pos,1),size(data_source.trial{iTrial},2));
     mapped_trial(leadfield.inside,:)=data_source.trial{iTrial};
     mapped_trial(~leadfield.inside,:)=NaN;
     mapped_trial(:,~leadfield.inside)=NaN;
     data_source.trial{iTrial}=mapped_trial;
    end
   end
   
   % get back trial markings
   if exist('trial_markings_sampleinfo','var')
    data_source.trial_markings_sampleinfo=trial_markings_sampleinfo;
   end
   if exist('trial_markings','var')
    data_source.trial_markings=trial_markings;
   end
   
   % save data
   save(result_file,'data_source','subject')
   
  else
   fprintf('LCMV: Broadband data present\n')
  end         
 otherwise
  error('No legitimate beamformer method')
end

  
%% if we got this far, place a stamp for completed processing
subject=nmri_stamp_subject(subject,['project_timecourse_' params.headmodel '_' params.src_project.beamformer] ,params);

end % final end

