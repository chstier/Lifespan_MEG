function nmri_write_mgh( fname, vox2ras ,metrics )
%nmri_write_mgh( fname, vox2ras ,metrics )
%   saves metrics as MGH for Freesurfer / PALM

if iscell(metrics)
 % make to double
 sz=size(metrics{1});
 if length(sz)>1 && sz(1)<sz(2)
  % transpose to have 1st dim max
  metrics=cellfun(@transpose,metrics,'un',0);
  sz=size(metrics{1});
 end
 % make 3 dims
 while length(sz)<3
  sz=[sz 1];
 end
 ml=length(metrics);
 sz=[sz ml];
 fmetrics=zeros(sz);
 for i=1:ml
  % now parse all
  fmetrics(:,:,:,i)=metrics{i};
 end
end
if isnumeric(metrics) 
 if length(size(metrics))~=4
  error('Metric seems to be numeric array, but does not have 4 dimensions. This is not supported. Give cell array for flexibility.')
 else
  fmetrics=metrics;
  clear metrics
 end
end


% now do the actual safe
save_mgh(fmetrics, fname, vox2ras);

end

