function [ out ] = nf_update_struct( orig, update )
%[ out ] = nf_update_struct( orig, update )
%   will update all fields in orig by respective fields in update
%   

items=fieldnames(update);
for i=1:length(items)
 if isfield(orig,items{i})
  if isstruct(orig.(items{i}))&&isstruct(update.(items{i}))
   orig.(items{i})=nf_update_struct(orig.(items{i}),update.(items{i}));
  else
   orig.(items{i})=update.(items{i});
  end
 else
  orig.(items{i})=update.(items{i});
 end
end
out=orig;
end

