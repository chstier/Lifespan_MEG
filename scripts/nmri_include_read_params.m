
% this will deal with the params structs only, subject is not used
% intended as matlab include


% check the params struct
if (~exist('params','var') || isempty(params))
 % check if the params are attached in the subject struct (this can be
 % needed in compiled mode)
 if (~exist('analysis_params.m','file'))
  error('Need to find analysis paramter file (analysis_params.m) in the current path, or have it in the call ')
 else
  analysis_params
  if (~exist('params','var')) 
   error('Problems with loading the parameter file (analysis_params.m)')  
  end
 end
end


    
