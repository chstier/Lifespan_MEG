function [ params ] = nmri_get_modality_params( params, mod_type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% check if we have a specific subfield
if (isfield(params,mod_type))

 params_b=params.(mod_type);
 params_a=rmfield(params,mod_type);

 
 % will use the versatible update_struct now
 params=nf_update_struct(params_a,params_b);
 
%  % then do the works
%  fn_a=fieldnames(params_a);
%  fn_b=fieldnames(params_b);
%  fn_conc=[fn_a;fn_b];
%  
%  [fn_all, fn_ind]=unique(fn_conc); 
%  % now re-sort to original
%  fn_all=fn_conc(sort(fn_ind));
%  vals=cell(length(fn_all),1); 
 
%  for i=1:length(fn_all)
%   % check if present in all
%   k_a=find(strcmp(fn_a,fn_all{i}));
%   k_b=find(strcmp(fn_b,fn_all{i}));
%   
%   if ~isempty(k_a) && ~isempty(k_b)
%    % present in both, concat if multiple members
%    vals{i}=catstruct(params_a.(fn_all{i}),params_b.(fn_all{i}));
%   elseif ~isempty(k_a) 
%    vals{i}=params_a.(fn_all{i});
%   elseif ~isempty(k_b) 
%    vals{i}=params_b.(fn_all{i});
%   else
%    error('We should never end here...')
%   end
%  end
%  params=cell2struct(vals,fn_all);
% end

end

