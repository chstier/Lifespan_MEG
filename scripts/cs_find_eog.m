function [subject] = cs_find_eog(subject)
% This function will load the muscle-cleaned data set (containing bad and good trial markings)
% and preprocessed EOG-data. After selecting only good trials in the MEG and EOG data, 
% EOG-ICA component will be selected based on amplitude and phase information.
% Inspired by mne-tools and Dammers et al. (2008; DOI 10.1109/TBME.2008.926677)
% Selection of ecg-components will be stored (eog_mag/grad_comp-file).

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

%% split dataset into MAG / GRAD and load EOG
% take only the good trial selection, which should be identical for MAGs and GRADs.
% call the subject and params include
nmri_include_read_ps
% call the central selection function
badTrials = [];
goodTrials = [1:length(data.trial)]; % start with all good
[goodTrials, badTrials] = nmri_trial_selector(subject,data,params);

cfg = [];
cfg.channel = {'MEG*1'}; % magnetometer
cfg.trials = goodTrials;
data_mag = ft_selectdata(cfg, data);

cfg = [];
cfg.channel = {'MEG*2', 'MEG*3'}; % gradiometer
cfg.trials = goodTrials;
data_grad = ft_selectdata(cfg, data);

% load physio-channels and select only good trials
if (~isfield(subject,'dws_filt_physio'))
 subject.dws_filt_physio = fullfile(subject.analysis_dir,subject.id,'processed',['dws_filt_physio_' subject.id '_' subject.exam_id '.mat']);
end

physio = load(subject.dws_filt_physio);

cfg = [];
cfg.channel    = {'EOG061'};
cfg.trials       = goodTrials;
eog1             = ft_selectdata(cfg, physio.data);
eog1.label{:}     = 'EOG1';

cfg              = [];
cfg.channel      = {'EOG062'};
cfg.trials       = goodTrials;
eog2             = ft_selectdata(cfg, physio.data);
eog2.label{:}     = 'EOG2';

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


%% bandpass filter eog and ICA components
cfg             = [];
cfg.bpfilter    = 'yes'; 
cfg.bpfreq      = [1 10];
cfg.bpfiltord   =  1;
cfg.bpfilttype  = 'but';
comp_mag.comp   = ft_preprocessing(cfg, comp_mag.comp);
comp_grad.comp   = ft_preprocessing(cfg, comp_grad.comp);
eog1       = ft_preprocessing(cfg, eog1);
eog2       = ft_preprocessing(cfg, eog2);

%% Amplitude-based: correlation over full data length
% redefine trials to have only 1 trial 
% needs trial structure which is given by the first and last sample of the
% data (trials during first cleaning are rejected already!2020/11/03).
% But to check, wether number of trials/samples are concordant in all MEG channels and eog, 
% this discards NA's in case there are any for missing data 

cfg     = [];
cfg.trl = [1 data_grad.sampleinfo(end,2) 0]; 
eog2_long = ft_redefinetrial(cfg,eog2);
eog1_long = ft_redefinetrial(cfg,eog1);
comp_mag_long.comp = ft_redefinetrial(cfg,comp_mag.comp);
comp_grad_long.comp = ft_redefinetrial(cfg,comp_grad.comp);

% take the sample values from comp.trial, transform to a matrix and get rid of NA's
all_dat_comp_mag = cell2mat(comp_mag_long.comp.trial);
all_dat_comp_grad = cell2mat(comp_grad_long.comp.trial);
all_dat_eog2 = cell2mat(eog2_long.trial);
all_dat_eog1 = cell2mat(eog1_long.trial);

comp_mag_clean= all_dat_comp_mag(:, all(~isnan(all_dat_comp_mag))); 
comp_grad_clean = all_dat_comp_grad(:, all(~isnan(all_dat_comp_grad)));
eog2_clean = all_dat_eog2(~isnan(all_dat_eog2));
eog1_clean = all_dat_eog1(~isnan(all_dat_eog1));

% compute correlation for MAG/GRAD components and eog1/eog2
mag_eog1 = {};
mag_eog2 = {};
grad_eog1 = {};
grad_eog2 = {};

for c = 1:length(comp_mag_long.comp.label)
 corr_mag_eog1 = abs(corrcoef(comp_mag_clean(c,:), eog1_clean)); 
 corr_mag_eog2 = abs(corrcoef(comp_mag_clean(c,:), eog2_clean));
 corr_grad_eog1 = abs(corrcoef(comp_grad_clean(c,:), eog1_clean)); 
 corr_grad_eog2 = abs(corrcoef(comp_grad_clean(c,:), eog2_clean));
 mag_eog1{c} = corr_mag_eog1(2);
 mag_eog2{c} = corr_mag_eog2(2);
 grad_eog1{c} = corr_grad_eog1(2);
 grad_eog2{c} = corr_grad_eog2(2);
end

corr_mag_eog1 = find(cell2mat(mag_eog1) > 0.4);
corr_mag_eog2 = find(cell2mat(mag_eog2) > 0.4);
corr_grad_eog1 = find(cell2mat(grad_eog1) > 0.4);
corr_grad_eog2 = find(cell2mat(grad_eog2) > 0.4);

% plot correlations 
% MAG components and eog1/eog2
corrFig1 = figure;
subplot(2,1,1); heatmap(cell2mat(mag_eog1),'XLabel','ICA components','YLabel','Trial (fulldata)','Title','correlation eog1-magnetometer');
subplot(2,1,2); heatmap(cell2mat(mag_eog2),'XLabel','ICA components','YLabel','Trial (fulldata)','Title','correlation eog2-magnetometer');
saveas(corrFig1,fullfile(subject.QCdir,['eog_mag_corr_' subject.id '_' subject.exam_id '.png']),'png');        
set(corrFig1,'Visible','on')
close

% GRAD components and eog1/eog2
corrFig2 = figure;
subplot(2,1,1); heatmap(cell2mat(grad_eog1),'XLabel','ICA components','YLabel','Trial (fulldata)','Title','correlation eog1-gradiometer');
subplot(2,1,2); heatmap(cell2mat(grad_eog2),'XLabel','ICA components','YLabel','Trial (fulldata)','Title','correlation eog2-gradiometer');
saveas(corrFig2,fullfile(subject.QCdir,['eog_grad_corr_' subject.id '_' subject.exam_id '.png']),'png');        
set(corrFig2,'Visible','on') 
close

clear mag_eog1 mag_eog2 grad_eog1 grad_eog2;


%% Now phase-based: coherence per trial
% append eog and ICA comp data
data_eog1_mag = ft_appenddata([],eog1, comp_mag.comp);
data_eog2_mag = ft_appenddata([],eog2, comp_mag.comp);
data_eog1_grad = ft_appenddata([],eog1, comp_grad.comp);
data_eog2_grad = ft_appenddata([],eog2, comp_grad.comp);

% compute a frequency decomposition of all components and the ECG
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.foilim     = [1 5];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq_mag1           = ft_freqanalysis(cfg, data_eog1_mag);
freq_mag2           = ft_freqanalysis(cfg, data_eog2_mag);
freq_grad1           = ft_freqanalysis(cfg, data_eog1_grad);
freq_grad2           = ft_freqanalysis(cfg, data_eog2_grad);

% compute coherence between all components and the ECG
cfg            = [];
cfg.channelcmb = {'all' 'EOG1'};
cfg.jackknife  = 'no';
cfg.method     = 'coh';
fdcomp_mag1         = ft_connectivityanalysis(cfg, freq_mag1);
fdcomp_grad1         = ft_connectivityanalysis(cfg, freq_grad1);

cfg            = [];
cfg.channelcmb = {'all' 'EOG2'};
cfg.jackknife  = 'no';
cfg.method     = 'coh';
fdcomp_mag2         = ft_connectivityanalysis(cfg, freq_mag2);
fdcomp_grad2        = ft_connectivityanalysis(cfg, freq_grad2);

% look at the coherence spectrum between all components and the ECG
% MAG
hFig1 = figure;
subplot(4,1,1); plot(fdcomp_mag1.freq, abs(fdcomp_mag1.cohspctrm)); legend('coherence eog1-magnetometer');
subplot(4,1,2); imagesc(abs(fdcomp_mag1.cohspctrm));
subplot(4,1,3); plot(fdcomp_mag2.freq, abs(fdcomp_mag2.cohspctrm)); legend('coherence eog2-magnetometer')
subplot(4,1,4); imagesc(abs(fdcomp_mag2.cohspctrm));

saveas(hFig1,fullfile(subject.QCdir,['eog_mag_coh_' subject.id '_' subject.exam_id '.png']),'png');        
set(hFig1,'Visible','on')
close

% GRAD
hFig2 = figure;
subplot(4,1,1); plot(fdcomp_grad1.freq, abs(fdcomp_grad1.cohspctrm)); legend('coherence eog1-gradiometer');
subplot(4,1,2); imagesc(abs(fdcomp_grad1.cohspctrm));
subplot(4,1,3); plot(fdcomp_grad2.freq, abs(fdcomp_grad2.cohspctrm)); legend('coherence eog2-gradiometer')
subplot(4,1,4); imagesc(abs(fdcomp_grad2.cohspctrm));

saveas(hFig2,fullfile(subject.QCdir,['eog_grad_coh_' subject.id '_' subject.exam_id '.png']),'png');        
set(hFig2,'Visible','on')
close

% compute the average coherence across frequencies
% MAG
av_coh_mag1 = mean(fdcomp_mag1.cohspctrm, 2);
coh_mag_eog1 = find(av_coh_mag1 > 0.3);

av_coh_mag2 = mean(fdcomp_mag2.cohspctrm, 2);
coh_mag_eog2 = find(av_coh_mag2 > 0.3);

% GRAD
av_coh_grad1 = mean(fdcomp_grad1.cohspctrm, 2);
coh_grad_eog1 = find(av_coh_grad1 > 0.3);

av_coh_grad2 = mean(fdcomp_grad2.cohspctrm, 2);
coh_grad_eog2 = find(av_coh_grad2 > 0.3);

%% Make final selection of components based on correlation and coherence 

mag_eog_comp = [corr_mag_eog1' corr_mag_eog2' coh_mag_eog1' coh_mag_eog2']; 
grad_eog_comp =  [corr_grad_eog1' corr_grad_eog2' coh_grad_eog1' coh_grad_eog2'];

% save MAG components
if (~isfield(subject,'mag_eog_comp'))
 subject.mag_eog_comp=fullfile(subject.analysis_dir,subject.id,'processed',['eog_mag_comp_' subject.id '_' subject.exam_id '.mat']);
end
save(subject.mag_eog_comp, 'mag_eog_comp')

% save GRAD components
if (~isfield(subject,'grad_eog_comp'))
 subject.grad_eog_comp=fullfile(subject.analysis_dir,subject.id,'processed',['eog_grad_comp_' subject.id '_' subject.exam_id '.mat']);
end
save(subject.grad_eog_comp, 'grad_eog_comp')

end
 %% plot ICA components
 % cfg          = [];
 % cfg.channel  = [4]; % components to be plotted
 % cfg.viewmode = 'component';
 % cfg.layout   = 'neuromag306mag_helmet.mat'; % specify the layout file that should be used for plotting
 % ft_databrowser(cfg, comp_mag.comp)
 % % 
 % cfg          = [];
 % cfg.channel  = [3]; % components to be plotted
 % cfg.viewmode = 'component';
 % cfg.layout   = 'neuromag306planar_helmet.mat'; % specify the layout file that should be used for plotting
 % ft_databrowser(cfg, comp_grad.comp)

% %% check whether subject has physio-data
% if length(data.label) == 309
% else
%  fprintf('Dataset does not contain ECG or EOG recordings. List subject.id and reject manually.')
%  file = 'physiodata_missing.mat';
%  if ~exist(fullfile(pwd, file))
%    missing_physio = {'ecg', 'eog'};
%    missing_physio{2,2} = subject.id;
%  else
%   load(fullfile(pwd, file));
%   if strcmp(missing_physio{length(missing_physio),1}, subject.id);
%    missing_physio{length(missing_physio),2} = subject.id;
%   else
%   missing_physio{length(missing_physio)+1,2} = subject.id;
%   end
%  end
%  save(file,'missing_physio')
% end