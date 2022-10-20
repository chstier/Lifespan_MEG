function [ params ] = nmri_get_params( analysis_dir )
%[ params ] = nmri_get_params( analysis_dir )
%   Returns the params for a specific analyis dir or pwd (if not specified)

if ~exist('analysis_dir','var')
 analysis_dir=pwd;
end
 
% check the params struct
if exist(analysis_dir,'dir')
 prev_dir=pwd;
 cd(analysis_dir)
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in analysis dir ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the paramter file (analysis_params.m)')  
  end
 end
 cd(prev_dir)
else
 error(['Analysis Dir ' analysis_dir ' not found'])
end

end

