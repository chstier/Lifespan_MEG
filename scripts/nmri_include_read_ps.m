
% this will deal with the subject and params structs and paths
% intended as matlab include



% now check the subject strcut
if ~isstruct(subject)
 % may be a mat file 
 if (exist(subject,'file'))
  if (strcmp(subject(end-1:end),'.m'))
   if isdeployed
    error('Cannot deal with a .m subject file in compiled/deployed exection. Provide a .mat file, or a struct.')
   else
    eval(subject);
   end
  elseif (strcmp(subject(end-1:end),'.mat'))
   load(subject);
  else
   error('Subject provided is not a valid struct, .m or .mat file')
  end
 else
  error('Subject provided is not a valid struct, .m or .mat file')
 end
else
 % we have a subject struct, but we should check for updates in the
 % subject_info.json file
 if (exist(fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id,'subject_info.json'),'file'))
  if ~exist('output','var') || ~islogical(output) || output
   disp('Have subject struct and found subject_info JSON file - loading and comparing');
  end
  osubject=subject;
  subject=jsondecode(fileread(fullfile(subject.analysis_dir,subject.id,'conf',subject.exam_id,'subject_info.json')));

  % now we need to check each field
  all_fields=fieldnames(subject);
  for i=1:length(all_fields)
   if (~isfield(osubject,all_fields{i}))
    % we do not know this one yet, so add
    osubject.(all_fields{i})=subject.(all_fields{i});
   else
    % Both have it, check equality
    if (~isequal(subject.(all_fields{i}),osubject.(all_fields{i})))
     if ~strcmp(all_fields{i},'analysis_dir')
      % unequal, overwrite unless analysis_dir
      if ~exist('output','var') || ~islogical(output) || output
       fprintf('Notice: overwriting subject struct (%s) from subject_info.json\n',all_fields{i})
      end
      osubject.(all_fields{i})=subject.(all_fields{i});
     end
    end
   end
  end
  % now write back
  subject=osubject;
  clear osubject;
 end
end

% check the params struct
if (~exist('params','var') || isempty(params))
 % check if the params are attached in the subject struct (this can be
 % needed in compiled mode)
 if isfield(subject,'params')
  params=subject.params;
 else
  if (~exist('analysis_params.m','file'))
   error('Need to find analysis paramter file (analysis_params.m) in the current path, or have it in the call ')
  else
   analysis_params
   if (~exist('params','var')) 
    error('Problems with loading the parameter file (analysis_params.m)')  
   end
  end
 end
end


% deal with layout
%if (~isfield(subject,'layout'))
% if (isfield(params,'layout'))
%  subject.layout=params.layout;
% elseif (isfield(params.(subject.dtype),'layout'))
%   subject.layout=params.(subject.dtype).layout;
% else  
%  error('no layout scheme specified in either subject or param')
% end
%end

    
