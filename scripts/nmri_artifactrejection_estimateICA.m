function [ subject, data ] = nmri_artifactrejection_estimateICA(subject, data, params)
%[ subject, data ] =  nmri_artifactrejection(subject, data, params)
%  
% This is a non-interactive function to estimate the ICA (gridable) 
%
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 09/2019



% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps


% check if we have dws_filt_dataset (this is the minimum)
if (~exist('data','var') ) 
 if (~isfield(subject,'clean_dataset') || ~exist(subject.clean_dataset,'file'))
  error('Visually cleaned dataset not found, run nmri_artifactrejection first')
 else
  load(subject.clean_dataset,'data');
 end
end

% make our output path and dir

% want ICA cleaning?
if (isfield(params,'useICA_clean') && params.useICA_clean==1)
 if (~isfield(subject,'cleanICA_dataset'))
  subject.cleanICA_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_' subject.id '_' subject.exam_id '.mat']);
 end
 if (~isfield(subject,'ICA_components'))
  subject.ICA_components=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_comp_' subject.id '_' subject.exam_id '.mat']);
 end
end


 %% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% now do the ICA cleaning if this is wanted
if (isfield(params,'useICA_clean') && params.useICA_clean==1)

 %% now do the ICA
 
 % temporarily reject technically bad trials
 rej_params=params;
 if isfield(rej_params,'rejectEvents')
  % do not reject events now, only stimuli and technical
  rej_params=rmfield(rej_params,'rejectEvents');
 end
 if isfield(rej_params,'rejectVigilance')
  % do not reject vigilance now, only stimuli and technical
  rej_params=rmfield(rej_params,'rejectVigilance');
 end
 
 [ goodTrials, ~ ] = nmri_trial_selector(subject,data,rej_params);
 
 data_r=data;
 tcfg           = [];
 tcfg.trials    = goodTrials;
 if isfield(data,'bad_channels')
 good_channels={};
 for i=1:length(data.label)
  if ~any(strcmp(data.label{i},data.bad_channels))
   good_channels(end+1)=data.label(i);
  end
 end
 tcfg.channel=good_channels;
end

 % remove trial markings to avoid Fieldtrip warnings
 if isfield(data_r,'trial_markings')
  data_r=removefields(data_r,'trial_markings');
 end
 if isfield(data_r,'trial_markings_sampleinfo')
  data_r=removefields(data_r,'trial_markings_sampleinfo');
 end
 if isfield(data_r,'bad_channels')
  data_r=removefields(data_r,'bad_channels');
 end
 
 data_r=ft_selectdata(tcfg,data_r);
  
  
  
 % set a safety margin of ICA componts...cannot be more than the number of
 % channels
 comp_offset=2;
  
 % use ICA in order to identify cardiac and blink components
 % this step can take a while
 cfg                 = [];
 cfg.method          = 'runica';
 cfg.runica.maxsteps = 128; % could be reduced for speed
 cfg.channel         = {'all','-refchan'};
 cfg.numcomponent    = min(100,(length(data_r.label)-comp_offset));
 w=warning('off','all');
 comp                = ft_componentanalysis(cfg, data_r);
 warning(w); % back to usual state  

 % check for complex numbers if ICA
 while (~isreal(comp.topo) && (comp_offset+10)<length(data.label))
  warning(['ICA failed for number of components=' num2str(cfg.numcomponent) '. Trying to lower by 2.'])
  comp_offset=comp_offset+2;
  % now try again with more offset
  cfg.numcomponent    = min(100,(length(data.label)-comp_offset));
  w=warning('off','all');
  comp                = ft_componentanalysis(cfg, data);
  warning(w); % back to usual state  
 end
 
 if ~isreal(comp.topo)
  error('ICA still failed...likely a problem with the data/number of channels')
 else
  % if we got this far, place a stamp for ICA estimation
  subject=nmri_stamp_subject(subject,'artifactrejectionICA_estimation',params);
  % save ICA data
  disp('Writing ICA data...')
  save(subject.ICA_components,'comp','subject')
 end

end



    

end

