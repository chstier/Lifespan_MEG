function [icv] = cs_get_icv(subject, suma_folder, mod)

%% This script loads the aseg.stats-file of each subject
% and extracts the estimated intracranial volume eTIV
% Christina Stier, written and modified in 2020/22

stats_dir = fullfile(suma_folder, [subject.id mod], 'stats');

% load aseg.stats-file 
fsstat_file = fullfile(stats_dir, 'aseg.stats');
 
fid = fopen(fullfile(stats_dir, 'aseg.stats'));
all_data = textscan(fid,'%s');
fclose(fid)

% take estimated intracranial volume eTIV
icv = all_data{1,1}{321,1};
icv = icv(1:end-8); % get rid of comma
icv = str2num(icv);
end