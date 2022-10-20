function [ filtered_subjects ] = nmri_filter_subjects ( all_subjects, filters )
%[ filtered_subjects ] = nmri_filter_subjects ( all_subjects, filters )
%   Filters subjects according to certain fields in the structs, very
%   flexible filters are possible
%
% all_subjects  = array of subject struct, e.g. from nmri_all_subjects
% filters       = struct array of filter
%
% if the struct field is logical, if will check for presense
% e.g. filters.stamps.processing_singleshell to check if this field is
% present / this processing step is done
%
% if the struct field is a char, it will match based on regular expression
% e.g. filters.exam_id='eeg/HD_EEG.rest.*' to match all HD_EEG.rest
% (see Matlab help, as in SPM)
%
% if the struct field is a cell array, each item will be checked and
% accepted if ANY of it matches. Usefull to filter a group of certain
% subjects.
% e.g. filters.id={'P.AABB01A','P.BBCC01.*'} would allow the subject
% P.AABB01A (timepoint A only) and all timepoints of P.BBCC01 
%               
% it is possible to recurse and combine. 
% e.g. filters.stamps.artifactrejection.user={'nfocke','jdoe'} would allow
% only subject being artifact rejected by nfocke or jdoe (maybe all the
% other raters were no good... ;) )
%
%

filtered_subjects=[];
for n=1:length(all_subjects)
 subject=all_subjects{n};
 if nmri_match_structs(subject, filters)
  filtered_subjects{end+1,1}=subject;
 end
end

end

