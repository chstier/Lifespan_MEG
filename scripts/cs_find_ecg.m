function [subject] = cs_find_ecg(subject)
% This function will load the muscle-cleaned data set (containing bad and good trial markings)
% and preprocessed ECG-data. After selecting only good trials in the MEG and ECG data, 
% coherence of ECG data and ICA components will be computed based on:
% http://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_ecg_artifacts/
% Selection of ecg-components will be stored (ecg_mag/grad_comp-file).

% Written and modified by Christina Stier, 2020/22

%% load cleaned dataset (all trials)
if (~isfield(subject,'clean_dataset'))
 subject.clean_dataset=fullfile(subject.analysis_dir,subject.id,'processed',['clean_' subject.id '_' subject.exam_id '.mat']);
end

load(subject.clean_dataset);

%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end

%% split dataset into MAG / GRAD and load ECG
% take only the good trial selection, which should be identical for MAGs and GRADs.
% call the subject and params include
nmri_include_read_ps
% call the central selection function
badTrials = [];
goodTrials = [1:length(data.trial)]; % start with all good
[goodTrials, badTrials] = nmri_trial_selector(subject,data,params);

cfg = [];
cfg.channel = {'MEG*1'};
cfg.trials = goodTrials;
data_mag = ft_selectdata(cfg, data);

cfg = [];
cfg.channel = {'MEG*2', 'MEG*3'};
cfg.trials = goodTrials;
data_grad = ft_selectdata(cfg, data);

% load physio-channels and select only good trials
if (~isfield(subject,'dws_filt_physio'))
 subject.dws_filt_physio = fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_physio_' subject.id '_' subject.exam_id '.mat']);
end

physio = load(subject.dws_filt_physio);

cfg              = [];
cfg.channel      = {'ECG*'};
cfg.trials       = goodTrials;
ecg              = ft_selectdata(cfg, physio.data);
ecg.label{:}     = 'ECG';

%% now load ICA components of MAGs and GRADs
if (~isfield(subject,'ICA_components_mag'))
  subject.ICA_components_mag=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_comp_mag_' subject.id '_' subject.exam_id '.mat']);
end

disp('Load ICA comp mag file...')
comp_mag = load(subject.ICA_components_mag, 'comp');

if (~isfield(subject,'ICA_components_grad'))
  subject.ICA_components_grad=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_comp_grad_' subject.id '_' subject.exam_id '.mat']);
end

disp('Load ICA comp grad file...')
comp_grad= load(subject.ICA_components_grad, 'comp');

%% detect QRS in ECG channel
cfg                       = [];
cfg.continuous            = 'yes';
cfg.artfctdef.ecg.pretim  = 0.25;
cfg.artfctdef.ecg.psttim  = 0.50-1/1500;
cfg.channel               = {'ECG'};
cfg.artfctdef.ecg.inspect = {'ECG'};
cfg.artfctdef.ecg.feedback = 'no';
[cfg, artifact]           = ft_artifact_ecg(cfg, ecg);

%%
% % bandpass filter 1-20Hz
% cfg             = [];
% cfg.bpfilter    = 'yes'; 
% cfg.bpfreq      = [1 20];
% cfg.bpfiltord   =  1;
% cfg.bpfilttype  = 'but'
% comp_filt   = ft_preprocessing(cfg, comp);
% ecg_filt        = ft_preprocessing(cfg, ecg);

%% define trials in ICA components according to ecg artifacts
cfg            = [];
cfg.trl        = [artifact zeros(size(artifact,1),1)];
comp_art_mag   = ft_redefinetrial(cfg, comp_mag.comp);
comp_art_grad  = ft_redefinetrial(cfg, comp_grad.comp);
ecg_art       = ft_redefinetrial(cfg, ecg);

data_all_mag = ft_appenddata([],ecg_art, comp_art_mag);
data_all_grad = ft_appenddata([],ecg_art, comp_art_grad);

% average the components timelocked to the QRS-complex
cfg           = [];
timelock_mag      = ft_timelockanalysis(cfg, data_all_mag);
timelock_grad      = ft_timelockanalysis(cfg, data_all_grad);

timelock = timelock_grad;

%% Plot timelocked data
hFig1 = figure;
subplot(3,1,1); plot(timelock_mag.time, timelock_mag.avg(1,:))
subplot(3,1,2); plot(timelock_mag.time, timelock_mag.avg(2:end,:))
subplot(3,1,3); imagesc(timelock_mag.avg(2:end,:));

saveas(hFig1,fullfile(subject.QCdir,['ecg_mag_timelock_' subject.id '_' subject.exam_id '.png']),'png');
set(hFig1,'Visible','on')

hFig2 = figure;
subplot(3,1,1); plot(timelock_grad.time, timelock_grad.avg(1,:))
subplot(3,1,2); plot(timelock_grad.time, timelock_grad.avg(2:end,:))
subplot(3,1,3); imagesc(timelock_grad.avg(2:end,:));

saveas(hFig2,fullfile(subject.QCdir,['ecg_grad_timelock_' subject.id '_' subject.exam_id '.png']),'png');
set(hFig2,'Visible','on')

%% compute coherence between ECG channel and ICA components
% compute a frequency decomposition of all components and the ECG
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.foilim     = [0 48];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq_mag           = ft_freqanalysis(cfg, data_all_mag);
freq_grad          = ft_freqanalysis(cfg, data_all_grad);

% compute coherence between all components and the ECG
cfg            = [];
cfg.channelcmb = {'all' 'ECG'};
cfg.jackknife  = 'no';
cfg.method     = 'coh';
fdcomp_mag     = ft_connectivityanalysis(cfg, freq_mag);
fdcomp_grad    = ft_connectivityanalysis(cfg, freq_grad);

%% look at the coherence spectrum between all components and the ECG
hFig3 = figure;
subplot(4,1,1); plot(fdcomp_mag.freq, abs(fdcomp_mag.cohspctrm)); legend('coherence ecg-magnetometer');
subplot(4,1,2); imagesc(abs(fdcomp_mag.cohspctrm));
subplot(4,1,3); plot(fdcomp_grad.freq, abs(fdcomp_grad.cohspctrm)); legend('coherence ecg-gradiometer')
subplot(4,1,4); imagesc(abs(fdcomp_grad.cohspctrm));

saveas(hFig3,fullfile(subject.QCdir,['ecg_mag_grad_coh_' subject.id '_' subject.exam_id '.png']),'png');        
set(hFig3,'Visible','on')

%% get the average coherence across frequencies
% magnetometer
av_coh_mag = mean(fdcomp_mag.cohspctrm, 2);
mag_ecg_comp = (find(av_coh_mag > 0.3))';

% gradiometer
av_coh_grad = mean(fdcomp_grad.cohspctrm, 2);
grad_ecg_comp = (find(av_coh_grad > 0.3))';

%% Get range for timelocked data, pick 2 channels with the highest range in case coherence did not work out

% for ICA components of the Magnetometer
if isempty(mag_ecg_comp) && isempty(grad_ecg_comp)
 ecg_file = fullfile(subject.analysis_dir, 'subjects_ecg_comp_missing.mat');
 
 if exist(ecg_file)  
  load(ecg_file);
  reason = 'coh did not work for ecg in both, mag and grad comp';
  ecg_subjects{end+1,1} = subject.id;
  ecg_subjects{end,2} = reason;
 else
  reason = 'coh did not work for ecg in both, mag and grad comp';
  ecg_subjects = {};
  ecg_subjects{end+1,1} = subject.id;
  ecg_subjects{end,2} = reason;
 end
 save('subjects_ecg_comp_missing.mat', 'ecg_subjects')
 
elseif isempty(mag_ecg_comp) || isempty(grad_ecg_comp)
 ecg_file = fullfile(subject.analysis_dir, 'subjects_ecg_comp_missing.mat'); 

 if exist(ecg_file)  
  load(ecg_file);
  reason = 'coh did not work for ecg in mag OR grad comp';
  ecg_subjects{end+1,1} = subject.id;
  ecg_subjects{end,2} = reason;
 else
  reason = 'coh did not work for ecg in mag OR grad comp';
  ecg_subjects = {};
  ecg_subjects{end+1,1} = subject.id;
  ecg_subjects{end,2} = reason;
 end
 save('subjects_ecg_comp_missing.mat', 'ecg_subjects')

end

if isempty(mag_ecg_comp)
% range_av_time_mag = {};
% for m = 1:length(timelock_mag.avg(:,1))s
%  range_av_time_mag{m,1} = range(timelock_mag.avg(m,:));
% end
% range_av_time_mag = cell2mat(range_av_time_mag(2:end)); % exclude the first row which is the ecg channel itself
% comp_ecgs = find(range_av_time_mag == max(range_av_time_mag));
peak_av_time_mag = {};
peak_av_time_mag = max(abs(timelock_mag.avg), [], 2); % get peak for each channel (absolute number as peaks can be negative)

peak_av_time_mag = peak_av_time_mag(2:end); % exclude ecg itself
peaks_sorted = sort(peak_av_time_mag); % sort channels based on peak (positive peaks here)

peak_comp_ecg1 = find(peak_av_time_mag == max(peak_av_time_mag)); % get the largest peak
peak_comp_ecg2 = find(peak_av_time_mag == peaks_sorted(end-1)); % get the 2nd largest peaks

mag_ecg_comp = [peak_comp_ecg1 peak_comp_ecg2];

clear peak_comp_ecg1 peak_comp_ecg2 peaks_sorted

end

% for ICA components of the Gradiometer

if isempty(grad_ecg_comp)
peak_av_time_grad = {};
peak_av_time_grad = max(abs(timelock_grad.avg), [], 2); % get peak for each channel (absolute number as peaks can be negative)

peak_av_time_grad = peak_av_time_grad(2:end); % exclude ecg itself 
peaks_sorted = sort(peak_av_time_grad); % sort channels based on peak

peak_comp_ecg1 = find(peak_av_time_grad == max(peak_av_time_grad)); % get the largest peak
peak_comp_ecg2 = find(peak_av_time_grad == peaks_sorted(end-1)); % get the 2nd largest peaks

grad_ecg_comp = [peak_comp_ecg1 peak_comp_ecg2];
end

%% save selections of ecg-components
if (~isfield(subject,'mag_ecg_comp'))
 subject.coh_selected=fullfile(subject.analysis_dir,subject.id,'processed',['ecg_mag_comp_' subject.id '_' subject.exam_id '.mat']);
end
save(subject.coh_selected, 'mag_ecg_comp')

if (~isfield(subject,'grad_ecg_comp'))
 subject.coh_selected=fullfile(subject.analysis_dir,subject.id,'processed',['ecg_grad_comp_' subject.id '_' subject.exam_id '.mat']);
end
save(subject.coh_selected, 'grad_ecg_comp')


close all
end

% %% plot ICA components
% cfg          = [];
% cfg.channel  = [6 12]; % components to be plotted
% cfg.viewmode = 'component';
% cfg.layout   = 'neuromag306mag_helmet.mat'; % specify the layout file that should be used for plotting
% ft_databrowser(cfg, comp_mag.comp)
% 
% cfg          = [];
% cfg.channel  = [22 27 36]; % components to be plotted
% cfg.viewmode = 'component';
% cfg.layout   = 'neuromag306planar_helmet.mat'; % specify the layout file that should be used for plotting
% ft_databrowser(cfg, comp_grad.comp)

 %% check whether subject has physio-data
% if length(data.label) == 309
% else
%   fprintf('Dataset does not contain ECG or EOG recordings. List subject.id and reject manually.')
%   file = 'physiodata_missing.mat';
%   if ~exist(fullfile(pwd, file))
%    missing_physio = {'ecg', 'eog'};
%    missing_physio{2,1} = subject.id;
%    save(file,'missing_physio')
%   else
%    load(fullfile(pwd, file));
%    missing_physio{length(missing_physio)+1,1} = subject.id;
%    end
%    save(file,'missing_physio')
%   end 
% end

%% check frequency decompositions for Grad / Mag
% cfg = [];
% cfg.output  = 'pow';
% cfg.channel = 'run*';
% cfg.method  = 'mtmfft';
% cfg.taper   = 'dpss';
% cfg.pad = 'nextpow2';
% cfg.foi     = 0.5:1:150; % 1/cfg1.length  = 1;
% cfg.tapsmofrq = 2;
% freq_grad_pow   = ft_freqanalysis(cfg, comp_grad.comp);
% 
% hFig = figure;
% hold on;
% plot(freq_grad_pow.freq, freq_grad_pow.powspctrm(:,:))
% legend('XX')
% xlabel('Frequency (Hz)');
% ylabel('Absolute Power (uV^2)');
