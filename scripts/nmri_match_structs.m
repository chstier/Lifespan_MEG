function [ match ] = nmri_match_structs( input, filters )
%function [ match ] = nmri_match_structs( input, filters )
%   will check if a struct (e.g. subject) matches with filter, allows to
%   recurse

items=fieldnames(filters);
match=true;
for i=1:length(items)
 if match && isfield(input,items{i}) 
  % we have it, so see what we want
  if ischar(filters.(items{i}))
   % char, regexpr match
   if ~ischar(input.(items{i})) || isempty(regexp(input.(items{i}),filters.(items{i})))
    match=false;
   end
  elseif isstruct(filters.(items{i}))
   if isstruct(input.(items{i})) 
    % both structs, so recurse
    match=nmri_match_structs(input.(items{i}),filters.(items{i}));
   else
    match=false;
   end
  elseif iscell(filters.(items{i}))
   % cell array, parse thorugh and check for any match
   any_match=false;
   for ii=1:length(filters.(items{i}))
    this_filter=[];
    this_filter.(items{i})=filters.(items{i}){ii};
    any_match=any_match||nmri_match_structs(input,this_filter);
   end
   if ~any_match
    match=false;
   end
  elseif islogical(filters.(items{i}))
   % just be happy with presence, already check above
  else
   % nothing fits, warn
   warning(sprintf('Filter item=%s was not parsed correctly. Should not happen, check filter definition',items{i}))
   match=false;
  end
 else
  match=false;
 end
end

end

