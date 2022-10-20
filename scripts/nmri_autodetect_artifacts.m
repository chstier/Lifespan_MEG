function [ subject, data ] = nmri_autodetect_artifacts( subject, data, params )
%[ subject, data ] = nmri_autodetect_artifacts( subject, data, params )
 
% This function will do automated trials and channel artifact selection
% similar to the vft_rejectvisual
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   data to analyze
% params    =   optional

% written by NF 09/2019 


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

% load data if needed
if (~exist('data','var') ) 
 % check if we have dws_filt_dataset (this is the minimum)
 if (~isfield(subject,'dws_filt_dataset') || ~exist(subject.dws_filt_dataset,'file'))
  error('Filtered and downsampled dataset not specified, run nmri_preproc first')
 else
  disp('Loading raw dataset...')
  load(subject.dws_filt_dataset,'data');
 end
end

%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end



%% make sure we have trial markings
if ~isfield(data,'trial_markings')
 data.trial_markings=cell(length(data.trial),4);
 % col1: sleep
 % col2: technical
 % col3: event
 % col4: rest/stimulation
 
 data.trial_markings_sampleinfo=cell(length(data.trial),2);
 % col1: sampleinfo start/stop
 % col2: seconds start/stop
 for i=1:length(data.trial)
  data.trial_markings_sampleinfo{i,1}=data.sampleinfo(i,:);
  data.trial_markings_sampleinfo{i,2}=data.sampleinfo(i,:)/data.fsample;
 end
end

if size(data.trial_markings,2)<4
 % make sure we have enough columns (backward compatabilty)
 for i=2:4
  if size(data.trial_markings,2)<i
   data.trial_markings=[data.trial_markings cell(size(data.trial_markings,1),1)];
  end
 end
end




%% now do the cleaning
bad_trials=zeros(1,length(data.trial));
bad_channels=zeros(length(data.label),1);

  
if isfield(params,'AUTO_clean') && isfield(params.AUTO_clean,'params')

 %% make a figure
 cols=0;
 for cn=1:length(params.AUTO_clean.params)
  cols=max([cols (length(params.AUTO_clean.params(cn).metrics)*2)]);
 end
 rows=length(params.AUTO_clean.params);
 hFig=figure('Position',[0,0,cols*300,rows*300],'Visible','on');
 
 % loop for all iterations
 for cn=1:length(params.AUTO_clean.params)
  
  if strcmp(params.AUTO_clean.params(cn).mode,'channels')
   opmode=2; % channels
   optxt='channels';
  elseif strcmp(params.AUTO_clean.params(cn).mode,'trials_by_channels')
   opmode=3; %trials-by-channels
   optxt='trials-by-channels';
  elseif strcmp(params.AUTO_clean.params(cn).mode,'trials')
   opmode=1; %trials
   optxt='trials';
  else
   error(['Unsupported mode=' params.AUTO_clean.params(cn).mode])
  end

  for n=1:length(params.AUTO_clean.params(cn).metrics)
   fprintf('Rejecting %s with metric=%s and cutoff=%d MAD\n',optxt,params.AUTO_clean.params(cn).metrics{n},params.AUTO_clean.params(cn).cutoff(n))
   level=calc_metric(data,params.AUTO_clean.params(cn).metrics{n},bad_trials,bad_channels);
  
 
   
   ocfg          = [];
ocfg.method   = 'summary'; 
if (isfield(subject,'layout'))
 ocfg.layout  = subject.layout ;
elseif (strcmp(subject.dtype,'MEG'))
 % preset for MEG
 ocfg.layout    = 'Neuromag306.lay';
else
 error('Layout could not be found...')
end



tcfg           = [];
tcfg.trials    = ~bad_trials;
tcfg.channel=data.label(~bad_channels);

data_r=data;

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


data_r         = ft_rejectvisual(ocfg,data_r); % press quit when done

clear data_r
  
   
   
   
   
   
   %bad=isoutlier(nanmean(level(:,:),opmode),'median','ThresholdFactor',params.AUTO_clean.params(cn).cutoff(n));
   
   % make our own
   if opmode==3
    %variance along trial axis
    vals=nanstd(level(:,:), [], 1).^2;
   else
    vals=nanmean(level(:,:),opmode);
   end
   MAD=mad(vals);
   MED=median(vals);
   cutoff=MED+(MAD*params.AUTO_clean.params(cn).cutoff(n));
   bad=vals>cutoff;
   
   fprintf('Median=%d, MAD=%d, Cutoff=%d\n',MED,MED,cutoff)
   
   
   
   subplot(rows,cols,(n*2)-1+(cols*(cn-1)))
   if opmode==2
    imagesc(level')
   else
    imagesc(level)
   end
   title([params.AUTO_clean.params(cn).metrics{n}])  
   
   subplot(rows,cols,(n*2)+(cols*(cn-1)))
   if opmode==1 || opmode==3
    hold on
    X=find(~bad_trials);
    plot(X(bad),vals(bad),'.','Color','r')
    plot(X(~bad),vals(~bad),'.','Color','b')
    line([1 numel(bad_trials)*1.1],[cutoff cutoff],'Color','r','LineStyle','--','LineWidth',1)
    xlabel('trials')
    hold off
   elseif opmode==2
    hold on
    X=find(~bad_channels);
    plot(X(bad),vals(bad),'.','Color','r')
    plot(X(~bad),vals(~bad),'.','Color','b')
    line([1 numel(bad_channels)*1.1],[cutoff cutoff],'Color','r','LineStyle','--','LineWidth',1)
    xlabel('channels')
    hold off
   end   
   title(['>' num2str(params.AUTO_clean.params(cn).cutoff(n))])  



   
   if opmode==1
    idx=find(~bad_trials); %these are the item checked
    bad_trials(idx)=bad_trials(idx)|bad;
    fprintf('...marking %d %s as bad, %d %s remaining\n',sum(bad),optxt,sum(~bad_trials),optxt)
    
    % debug
    for i=1:numel(bad)
     if bad(i)
      fprintf('Rejected trial=%d, val=%d, Cutoff=%d\n',idx(i),vals(i),cutoff)
     end
    end
    
   elseif opmode==3
    %variance along trial axis
    idx=find(~bad_trials); %these are the item checked
    bad_trials(idx)=bad_trials(idx)|bad;
    fprintf('...marking %d %s as bad, %d %s remaining\n',sum(bad),optxt,sum(~bad_trials),optxt)
    
    % debug
    for i=1:numel(bad)
     if bad(i)
      fprintf('Rejected trial=%d, val=%d, Cutoff=%d\n',idx(i),vals(i),cutoff)
     end
    end
   elseif opmode==2
    idx=find(~bad_channels); %these are the item checked
    bad_channels(idx)=bad_channels(idx)|bad;
    fprintf('...marking %d %s as bad, %d %s remaining\n',sum(bad),optxt,sum(~bad_channels),optxt)
    
    % debug
    for i=1:numel(bad)
     if bad(i)
      fprintf('Rejected channel=%s, val=%d, Cutoff=%d\n',data.label{idx(i)},vals(i),cutoff)
     end
    end
   end   
   %title(['AFTER ' params.AUTO_clean.params(cn).metrics{n} '>' num2str(params.AUTO_clean.params(cn).cutoff(n))])

  for i=1:numel(bad_channels)
   if bad_channels(i)==1
    fprintf('Bad Channel=%s\n',data.label{i})
   end
  end
  
  for i=1:numel(bad_trials)
   if bad_trials(i)==1
    fprintf('Bad Trial=%d\n',i)
   end
  end
  
  
  end
  fprintf('-----------------\n')
 end

end




saveas(hFig,fullfile(subject.QCdir,[ 'AutoReject_' subject.id '_' subject.exam_id '.png']),'png'); 
% only show after save
set(hFig,'Visible','on')


% now put the info into data and reject channels
subject.artifactrejection_AUTO.bad_trials=bad_trials;
subject.artifactrejection_AUTO.bad_channels=data.label(bad_channels==1);


% add to bad channels in data
if ~isfield(data,'bad_channels')
 data.bad_channels={};
end
for i=1:length(data.label)
 if bad_channels(i)
  fprintf('Rejecting Channel=%s\n',data.label{i})
  if ~any(strcmp(data.bad_channels,data.label{i}))
   data.bad_channels(end+1)=data.label(i);
  end
 end
end






% mark trials as bad
for i=1:length(data.trial)
 if bad_trials(i)
  % this trials seems to have been rejected - mark as bad (false)
  data.trial_markings{i,2}=false;  
 end
end




end


%% from ft_rejectvisual / rejectvisual_summary
% info
function level=calc_metric(data,metric,bad_trials,bad_channels)
 info=[];
 info.nchan=sum(~bad_channels);
 info.ntrl=sum(~bad_trials); 
 level = zeros(info.nchan, info.ntrl);
 idx=find(~bad_trials); % skip bad trials
  
 if strcmp(metric, 'zvalue') || strcmp(metric, 'maxzvalue') || strcmp(metric, 'abszvalue')
  % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(info.nchan, 1);
  runss  = zeros(info.nchan, 1);
  runnum = 0;

  for i=1:info.ntrl
    dat = data.trial{idx(i)}(~bad_channels,:); % skip bad channels
    runsum = runsum + nansum(dat, 2);
    runss  = runss  + nansum(dat.^2, 2);
    runnum = runnum + sum(isfinite(dat), 2);
  end
  mval = runsum./runnum;
  sd   = sqrt(runss./runnum - (runsum./runnum).^2);
 end

for i=1:info.ntrl
  dat = data.trial{idx(i)}(~bad_channels,:); % skip bad channels
   
  switch metric
    case 'var'
      level(:, i) = nanstd(dat, [], 2).^2;
    case 'min'
      level(:, i) = nanmin(dat, [], 2);
    case 'max'
      level(:, i) = nanmax(dat, [], 2);
    case 'maxabs'
      level(:, i) = nanmax(abs(dat), [], 2);
    case 'range'
      level(:, i) = nanmax(dat, [], 2) - nanmin(dat, [], 2);
    case 'kurtosis'
      level(:, i) = kurtosis(dat, [], 2);
    case '1/kurtosis'
      level(:, i) = 1./(kurtosis(dat, [], 2)+0.0000000000000000001); % make sure it is never 0
    case '1/var'
      level(:, i) = 1./((nanstd(dat, [], 2).^2)+0.0000000000000000001); % make sure it is never 0
    case 'zvalue'
      level(:, i) = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
    case 'maxzvalue'
      level(:, i) = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
    case 'abszvalue'
      level(:, i) = abs(nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2));
    otherwise
      ft_error('unsupported method');
  end
end


end


