function [ subject, data ] = cs_nmri_artifactrejection_reviewICA_mag(subject, data, params)
%[ subject, data ] =  nmri_artifactrejection(subject, data, params)
%  
% This is an interactive function to review and reject the ICA components 
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
 if (~isfield(subject,'clean_rej_dataset') || ~exist(subject.clean_rej_dataset,'file'))
  error('Visually cleaned dataset not found, run nmri_artifactrejection first')
 else
  load(subject.clean_rej_dataset,'data');
 end
end

% load the components
if (~isfield(subject,'ICA_components_mag') || ~exist(subject.ICA_components_mag,'file'))
 error('ICA components not found, run nmri_artifactrejection_estimateICA first')
else
 % check if has selected
 content=who('-file',subject.ICA_components_mag);
 if ismember('selected',content)
  load(subject.ICA_components_mag,'comp','selected');
 else
  load(subject.ICA_components_mag,'comp');
 end
end

if (isfield(params,'useICA_clean') && params.useICA_clean==1)
 if (~isfield(subject,'cleanICA_mag_dataset'))
  subject.cleanICA_mag_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['cleanICA_mag_' subject.id '_' subject.exam_id '.mat']);
 end
end
% 
% 
% % check for previous version
% if (isfield(params,'useICA_clean') && params.useICA_clean==1 && isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
%  % have ICA
%  button=questdlg({'A previously generated ICA-cleaned dataset was found.','This will be overwritten if you continue with the process.','ICA rejection always starts from the cleaned dataset.',...
%   'I.e. all edits of the previous file may be lost (e.g. later channel rejections)','Do you want to continue?'},'No');
%  if ~strcmpi(button,'Yes')
%   disp('Stopping on user request...')
%   return
%  end
% end


%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end




 %% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% now do the ICA cleaning if this is wanted
if (isfield(params,'useICA_clean') && params.useICA_clean==1)

 
 % visualize components in order to selected the components to be rejected
 cfg           = [];
 if (isfield(subject,'layout'))
  cfg.layout  = subject.layout ;
 elseif (strcmp(subject.dtype,'MEG'))
  % preset for MEG
  cfg.layout    = 'Neuromag360.lay';
 else
  error('Layout could not be found...')
 end
 
 if exist('selected','var')
  % set the previous selection
  cfg.selected=selected;
 end


 % use our own browser for component rejection
 h=msgbox({'Please select components to be rejected / kept'},'Visual ICA rejection');
 selected=ica_selection(cfg,comp);
 try
  % close, if not already
  delete(h)
 catch
 end
  
 % save ICA data
 save(subject.ICA_components_mag,'comp','selected')
         
 % remove the bad components and backproject the data
 cfg = [];
 cfg.component = find(~selected);
 if ~isempty(cfg.component)
 % w=warning('off','all');
 
 
  
  
  % need to deal with trials that are in data, but not in data_r/comp
  % because these were bad
  comp_r=comp;
  comp_r.trial=cell(size(data.trial));
  comp_r.time=cell(size(data.time));
  comp_r.sampleinfo=data.sampleinfo;
  
  
  for i=1:length(comp_r.trial)
   found_t=find(comp_r.sampleinfo(i,1)==comp.sampleinfo(:,1));
   if isempty(found_t)
    % this trials seems to have been rejected, put zeros for this
    % component
    comp_r.trial{i}=zeros(size(data.trial{i}));
    comp_r.time{i}=zeros(size(data.time{i}));
   else
    comp_r.trial{i}=comp.trial{found_t};
    comp_r.time{i}=comp.time{found_t};
   end
  end
  
  % Note: it seems that Fieldtrip is smart enough to deal with differences
  % in channels between comp and data
  
  data = ft_rejectcomponent(cfg, comp, data);
  %warning(w); % back to usual state
 end
  
 %% if we got this far, place a stamp for completed artifact rejection
 subject=nmri_stamp_subject(subject,'artifactrejectionICA_mag',params);
  
 %% save ICA cleaned data - after ICA
 nmri_write_dataset(subject.cleanICA_mag_dataset,data,subject);
end 


%% check if we want a QC power plot
if (isfield(params,'QC_plots') && params.QC_plots==1)
 % make plot
 
 button=questdlg({'Do you want to generate a QC power plot now?','Can take a few minutes'},'Generate QC plot?','No');
 if strcmpi(button,'Yes')
  
  % call the central QC power plot script
  cs_nmri_qcplot_power_mag(subject,data,params,['pow_mag_plot_postICA_' datestr(now,'yyyymmdd')]);
  
 end
end



    

end

