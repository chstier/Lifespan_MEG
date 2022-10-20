function [ subject ] = nmri_reroot_subject( subject, new_analysis_dir )
% [ subject ] = nmri_reroot_subject( subject )
%  Will fix a subject struct that was moved to a new location / analysis
%  dir
%
%  subject          = subject struct (input=old, output=re-rooted)
%  new_analysis_dir = you can provide the new analyis dir, if not set, pwd
%  is used

if ~exist('new_analysis_dir','var')
 new_analysis_dir=pwd;
end

if isfield(subject,'analysis_dir')
 old_analysis_dir=subject.analysis_dir;
 items=fieldnames(subject);
 for i=1:length(items) 
  if ischar(subject.(items{i}))
   subject.(items{i})=strrep(subject.(items{i}),old_analysis_dir,new_analysis_dir);
  end
 end
else
 warning('Could not read analysis dir from subject struct -- Probably incorrect struct provided')
end

