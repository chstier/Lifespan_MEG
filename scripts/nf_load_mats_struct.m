function [ concat ] = nf_load_mats_struct( mats, target_var )
%[ concat ] = nf_load_mats_struct( mats, target_var)
%   Load a cell array of mat files into a struct array
%   will fill up missing fields with []

if ~iscell(mats) || isempty(mats)
    error('Need call array with at least one file')
end
concat=[];
for i=1:length(mats)
    if ~exist(mats{i},'file')
        warning(['Mat file=' mats{i} ' not found'])
    else
     this_mat=load(mats{i},target_var);
     % add loaded_from stamp
     for fn=1:size(this_mat.(target_var),1)
        this_mat.(target_var)(fn).loaded_from=mats{i};
     end
     if isfield(this_mat,target_var)
         if ~isempty(concat)
             % concat all fields
             a=fieldnames(this_mat.(target_var));
             b=fieldnames(concat(1));
             ab=unique([a;b]);
             for f=1:length(ab)
                 if ~isfield(this_mat.(target_var),ab{f})
                    this_mat.(target_var)(1).(ab{f})=[];
                 end
                 if ~isfield(concat(1),ab{f})
                    concat(1).(ab{f})=[];
                 end
             end
             concat=[concat;this_mat.(target_var)];
         else
            concat=this_mat.(target_var);
         end
     end        
    end
end


end

