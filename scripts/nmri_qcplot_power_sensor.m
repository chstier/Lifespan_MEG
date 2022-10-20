function [ subject ] = nmri_qcplot_power_sensor(subject, data, params, title_txt)
%[ subject ] = nmri_qcplot_power_sensor(subject, data, params, title_txt)
%  
% Generates a QC plot for power at sensor space
%
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   data to be plotted
% title_txt =   if special title is asked for


% written by NF 09/2019



% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps


% check if we cleaned data
if (~exist('data','var') || isempty(data)) 
 if (isfield(subject,'cleanICA_dataset') && exist(subject.cleanICA_dataset,'file'))
  load(subject.cleanICA_dataset,'data');
 elseif (isfield(subject,'clean_dataset') && exist(subject.clean_dataset,'file'))
  load(subject.clean_dataset,'data');
 else 
  error('No data given and no cleaned dataset loadable')
 end
end

if (~exist('title_txt','var') || isempty(title_txt)) 
 title_txt='artifact_pow_sens_plot';
end


%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end



 %% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );



%% check if we want a QC power plot
if (isfield(params,'QC_plots') && params.QC_plots==1)
 % make plot
 

  % call the central selection function now
  [ goodTrials, badTrials ] = nmri_trial_selector(subject,data,params);

  fprintf('\nTotal: GoodTrials, N=%d / BadTrials, N=%d\n',length(goodTrials),length(badTrials))

  tcfg          = [];
  tcfg.trials   = goodTrials;    
  if isfield(data,'bad_channels')
   good_channels={};
   for i=1:length(data.label)
    if ~any(strcmp(data.label{i},data.bad_channels))
     good_channels(end+1)=data.label(i);
    end
   end
   tcfg.channel=good_channels;
   fprintf('\nTotal: GoodChannels, N=%d / BadChannels, N=%d\n',length(good_channels),length(data.bad_channels))
  else
   disp('No bad channels found')
  end



  % remove trial markings to avoid Fieldtrip warnings

  if isfield(data,'trial_markings')
   data=removefields(data,'trial_markings');
  end
  if isfield(data,'trial_markings_sampleinfo')
   data=removefields(data,'trial_markings_sampleinfo');
  end
  if isfield(data,'bad_channels')
   data=removefields(data,'bad_channels');
  end

  data=ft_selectdata(tcfg,data);

  for ff = 1:length(params.freqs)   
   fprintf('\nCalculation of frequency power for band=%s\n',params.freqsNames{ff})
   cfg           = [];
   cfg.method    = 'mtmfft';
   cfg.output    = 'pow';
   cfg.pad       = 'nextpow2'; % recommended for speed
   cfg.foi       = params.freqs(ff);
   cfg.tapsmofrq = params.tapsmofrq(ff); 
   pow_freq(ff)       = ft_freqanalysis(cfg, data);
  end

  % write the QC image

  disp('Now doing sensor space power plots for QC')   
  allf=length(params.freqs);
  cols=3;
  rows=ceil(allf/cols);
  hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','off'); 
  if (~isfield(subject,'layout'))
   if (isfield(params,'layout'))
    subject.layout=params.layout;
   else
    error('no layout scheme specified in either subject or param')
   end
  end
  % get our layout
  cfg = [];
  cfg.layout=subject.layout;
  cfg.skipscale='yes';
  cfg.skipcomnt='yes';
  layout=ft_prepare_layout(cfg);
  layout.label=upper(layout.label); % make upper case

  cfg=[];
  cfg.gridscale = 200; %more beautiful
  cfg.layout    = layout;
  cfg.comment   = 'xlim';
  for i=1:allf
  subplot(rows,cols,i)
  title(params.freqsNames{i},'FontSize',12,'FontWeight','bold')
  ft_topoplotTFR(cfg,pow_freq(i))
  end
  
  saveas(hFig,fullfile(subject.QCdir,[title_txt '_' subject.id '_' subject.exam_id '.png']),'png'); 
  % only show after save
  set(hFig,'Visible','on')


  %system(['eog ' fullfile(subject.QCdir,['artifact_pow_sens_plot_' subject.id '_' subject.exam_id '.png'])])

end



    

end

