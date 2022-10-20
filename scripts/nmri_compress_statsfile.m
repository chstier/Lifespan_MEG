function [ all_subjects ] = nmri_compress_statsfile( all_subjects , params )
%[ all_subjects ] = nmri_compress_statsfile( all_subjects )
%   Will shrink the stats-file to save disk space
%   will just retrain the power estmates

% check the params struct
if (~exist('params','var'))
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in the current path, or have it in the call ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the paramter file (analysis_params.m)')  
  end
 end
end

if isstruct(all_subjects)
 % in case a single subject is given
 all_subjects={all_subjects};
end


for n=1:length(all_subjects)
 % parse over all subjects
 subject=all_subjects{n};
 fprintf('Shrinking (%d/%d) Subject=%s, Exam_ID=%s ',n,length(all_subjects),subject.id,subject.exam_id)
 
 [ this_params ] = nmri_get_modality_params( params, subject.dtype );
 if exist(subject.stats,'file')
  % check size
  fstats = dir(subject.stats);
  if fstats.bytes>(10*1024*10124)
   % at least 10 mB
   load(subject.stats);
   source_filters_red=cell(1,length(source_filters));
   for i=1:length(source_filters)
    source_filters_red{i}.avg=source_filters{i}.avg;
    source_filters_red{i}.avg=removefields(source_filters_red{i}.avg,{'ori','eta','filter','label'});
   end
   source_filters=source_filters_red;
   % now overwrite
   save(subject.stats,'source_filters')
   fprintf('...shrunken\n')
  else
   fprintf('...already small - skipping\n\n')
  end
 else
  frprintf('Stat file=%s not found\n',subject.stats)
 end
end


end

