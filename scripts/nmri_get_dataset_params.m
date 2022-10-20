function [ params ] = nmri_get_dataset_params( params, dataset_mapping )
% will modify the params according to dataset mapping
if ~isempty(dataset_mapping)
 if ~isstruct(dataset_mapping)
  % should be a file then
  if ~exist(dataset_mapping,'file')
   error(['could not find ' dataset_mapping])
  end
  try
   ds_params=load(dataset_mapping,'params');
  catch
   warning(['could not find a params variable in ' dataset_mapping '\nWill assume this is intentionally'])
   ds_params=[];
  end

  if ~isempty(ds_params) && isfield(ds_params,'params')
   params=nf_update_struct(params,ds_params.params);
  end
 end
end

end

