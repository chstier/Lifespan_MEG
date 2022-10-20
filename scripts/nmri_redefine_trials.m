function [ subject, data ] = nmri_redefine_trials( cfg, subject, data )
%[ subject, data ] = nmri_redefine_trials( subject, data )
%   Will redefine trials (using fieldtrip) and take care of markings /
%   vigilance scoring, etc...
%   
% 06/2018 by NF



data2=ft_redefinetrial(cfg,data);

% now rebuild the trialmarkings
data2.trial_markings=cell(length(data2.trial),3);
% col1: sleep
% col2: technical
% col3: event

data2.trial_markings_sampleinfo=cell(length(data2.trial),2);
% col1: sampleinfo start/stop
% col2: seconds start/stop
for i=1:length(data2.trial)
 data2.trial_markings_sampleinfo{i,1}=data2.sampleinfo(i,:);
 data2.trial_markings_sampleinfo{i,2}=data2.sampleinfo(i,:)/data2.fsample;
end

% now fill the info
for i=1:length(data2.trial)
 % check which old trial this refers to - use seconds just in case
 startstop=data2.trial_markings_sampleinfo{i,2};
 foundtrial=[];
 for ii=1:length(data.trial) 
  if startstop(1)>=data.trial_markings_sampleinfo{ii,2}(1) && startstop(2)<=data.trial_markings_sampleinfo{ii,2}(2)
   foundtrial=ii;
   break
  end
 end
 data2.trial_markings(i,:)=data.trial_markings(foundtrial,:);
end

% update subject
subject.data_info.trial_markings=data2.trial_markings;
subject.data_info.trial_markings_sampleinfo=data2.trial_markings_sampleinfo;
subject.data_info.ntrials=length(data2.trial);

data=data2;

end


