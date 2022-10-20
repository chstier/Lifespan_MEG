function[subject,data] = cs_freqplot(subject,data)

% This fuction loads the ICA-cleaned dataset and performs frequency analysis 
% on GRADs and MAGs and saves frequency-QC plots (1-150Hz)
%
% Written and modified by Christina Stier 2020/2022

%% Load dws-filt_dataset
if (~exist('data','var') ) 
 % check if we have clean_rej_dataset (this is the minimum)
 if (~isfield(subject,'cleanICA_dataset') || ~exist(subject.cleanICA_dataset,'file'))
  error('Dataset not specified, run automatic artefactrejection (muslce) first')
 else
  disp('Loading cleanICA_dataset...')
  load(subject.cleanICA_dataset,'data');
 end
end

%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end

%% do fft-plot MEGGRAD
cfg = [];
cfg.output  = 'pow';
cfg.channel = 'meggrad';
cfg.method  = 'mtmfft';
cfg.taper   = 'dpss';
cfg.pad = 'nextpow2';
cfg.foi     = 0.5:1:150; % 1/cfg1.length  = 1;
cfg.tapsmofrq = 2;
data_freq   = ft_freqanalysis(cfg, data);

hFig = figure;
hold on;
plot(data_freq.freq, data_freq.powspctrm(:,:))
legend('All Gradiometer Neuromag306')
xlabel('Frequency (Hz)');
ylabel('Absolute Power (uV^2)');

%% save freq plot GRAD
if (~exist(fullfile(subject.QCdir,['grad_sensor_freqplot_' subject.id '_' subject.exam_id '.png']),'file'))
  disp('Now doing frequency plots for QC') 

  saveas(hFig,fullfile(subject.QCdir,['grad_sensor_freqplot_' subject.id '_' subject.exam_id '.png']),'png');        
  set(hFig,'Visible','on')
end

clear data_freq hFig

%% do fft-plot MEGMAG
cfg = [];
cfg.output  = 'pow';
cfg.channel = 'megmag';
cfg.method  = 'mtmfft';
cfg.taper   = 'dpss';
cfg.pad = 'nextpow2';
cfg.foi     = 0.5:1:150; % 1/cfg1.length  = 1;
cfg.tapsmofrq = 2;
data_freq   = ft_freqanalysis(cfg, data);

hFig = figure;
hold on;
plot(data_freq.freq, data_freq.powspctrm(:,:))
legend('All Magnetometer Neuromag306')
xlabel('Frequency (Hz)');
ylabel('Absolute Power (uV^2)');

%% save freq plot MEGMAG
if (~exist(fullfile(subject.QCdir,['mag_sensor_freqplot_' subject.id '_' subject.exam_id '.png']),'file'))
  disp('Now doing frequency plots for QC') 

  saveas(hFig,fullfile(subject.QCdir,['mag_sensor_freqplot_' subject.id '_' subject.exam_id '.png']),'png');        
  set(hFig,'Visible','on')
end

clear data_freq hFig
close all


