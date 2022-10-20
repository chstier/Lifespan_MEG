function [ subject ] = nmri_regenerate_cleaning(subject_new, subject_clean, params)
%[ subject ] = nmri_regenerate_cleaning(subject_new, subject_clean, params)
%  
% This function will re-generate the cleaning, marking and vigilance
% scoring from one subject struct to another subject (or the same one with
% updates name)
% Will use the dws-filt of subject_new and generate (new) cleaned datasets
% Note: ICA-cleaning will have to be re-run 
% 
% subject_new   =   subject structure of new subject, i.e. the one that was
%                   recently read, dws-filt
% subject_clean =   subject structure of cleaned subject, i.e. the one that was
%                   processed before. Can be the same as subject_clean, in
%                   which case the original files are overwritten
%
% params        =   anaylsis parameter struct (optional, will search for
%                    analysis_params.m if not set)

% written by NF 10/2019


% check the call
if (~exist('subject_new','var') ) 
 error('Need a valid subject struct to work with')
end

if (~exist('subject_clean','var') ) 
 fprintf('Only one subject provided\nCAVE: Overwriting data for subject=%s, exam=%s\n',subject.id,subject.exam_id)
 subject_clean=subject_new;
end





end

