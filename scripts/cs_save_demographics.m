%% This script extracts the demographics of selected subjects and saves 
% in a table in the following order: subj-id, age, sex, handedness
% Christina Stier, 2020/2022

% get selection of subjects

load('subjectnames_final.mat')

% get full list of all subjects available

subjlist = importdata('participants.tsv');

subjects = subjlist.textdata(2:end,1); % get subj names
subjlist.data(subjlist.data == 1) = zeros; % get index sex and change to 0 (males) and 1 (females)
subjlist.data(subjlist.data == 2) = ones;

subjects(:,2) = subjlist.textdata(2:end,2);
subjects(:,3) = num2cell(subjlist.data);
subjects(:,4) = subjlist.textdata(2:end,3);

% find demographics for those subjects needed

sel2 = ismember(subjects(:,1), subjlist_final);
final_sel = subjects(sel2, :);

% save final selection including demographics in a table
T = table(final_sel);
writetable(T, 'Demographics_cohortX_final.xlsx','WriteVariableNames',false)