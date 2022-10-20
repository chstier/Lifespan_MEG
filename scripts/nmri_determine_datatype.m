function [ subject ] = nmri_determine_datatype( subject )
%[ subject ] = nmri_determine_datatype( subject )
%   central script to determine datatype (EEG/MEG/...)
%   layout (for Fieldtrip viz and ICA)
%   montages (for eeg_score viewer)
%   will also try to deal with invasive EEG data


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct to work with')
end

% read hdr, if not there
if (~isfield(subject,'hdr') || isempty(subject.hdr))
 if strcmpi(subject.raw_dataset(end-3:end),'.mff')
  % check which MFF reader to use, clipped files are not read by v1...
  subject.hdr=ft_read_header(subject.raw_dataset,'headerformat',nmri_check_mff_reader(subject.raw_dataset));
  % need to reomve java Objects...these do not serialize
  if isfield(subject,'hdr') && isfield(subject.hdr,'orig') && isfield(subject.hdr.orig,'javaObjs')   
   subject.hdr.orig=rmfield(subject.hdr.orig,'javaObjs');
  end
  % also get rid of MFF_v3 original data...stupid idea to store that here
  if isfield(subject,'hdr') && isfield(subject.hdr,'orig') && isfield(subject.hdr.orig,'data')   
   subject.hdr.orig=rmfield(subject.hdr.orig,'data');
  end
  
  
 else
  subject.hdr=ft_read_header(subject.raw_dataset);
 end
end

% remove lfp label
lfp=strcmpi(subject.hdr.chantype,'lfp');
if sum(lfp)>1
 subject.hdr.chantype(lfp)=repmat({'unknown'},[sum(lfp) 1]);
end

% may need to deal with chantype of invasice EEG/ECoG that are mapped as
% unknown
if ((sum(strcmpi(subject.hdr.chantype,'eeg')) + sum(~cellfun(@isempty,regexpi(subject.hdr.chantype,'meg')))) / sum(strcmpi(subject.hdr.chantype,'unknown'))) < 0.5
 % at least 50% unknown channels, may be invasive data
 disp('Found many unknown channels, will assume this is invasive data')
 unk=subject.hdr.label(strcmpi(subject.hdr.chantype,'unknown'));
 repN=0;
 for i=1:length(unk)
  % check if still unknown
  if strcmpi(subject.hdr.chantype{strcmp(subject.hdr.label,unk{i})},'unknown')
   alpha=isstrprop(unk{i},'alpha');
   if alpha(1) && ~strcmp(unk{i}(alpha),'DC') %starts with a letter and not DC
    %find possible partners
    m=~cellfun(@isempty,regexp(unk,['^' unk{i}(alpha) '[0-9]*$']));
    if sum(m)>4 % smallest elec is 5 contacts for DIXI
     % think this is EEG and microVolts
     subject.hdr.chantype(m)=repmat({'eeg'},[sum(m) 1]);
     subject.hdr.chanunit(m)=repmat({'uV'},[sum(m) 1]);
     repN=repN+sum(m);
    end   
   end
  end
 end
 if repN>10
  subject.dtype='EEG_invasive';
  subject.repN=repN;
 end
end

%% now read our dataset mappings
mapping_dir=fullfile(subject.analysis_dir,'conf','dataset-mappings');
if ~exist(mapping_dir,'dir')
 error(['Could not find the dataset-mappings in ' mapping_dir])
end
all_d=dir(mapping_dir);
all_m={};
for i=1:size(all_d,1)
 if ~strcmp(all_d(i).name(1),'.') && ~all_d(i).isdir
  all_m=[all_m {fullfile(mapping_dir,all_d(i).name)}];
 end
end
dataset_id=nf_load_mats_struct(all_m,'dataset_id');

% now determine, what we have and set subject fields corresponding
best_match=0;
best_i=[];
% loop for all presets
for i=1:size(dataset_id,1)
 match=0;
 % loop for channels
 for ch=1:length(subject.hdr.label)
  if any(strcmpi(dataset_id(i).label,subject.hdr.label{ch}))
   match=match+1;
  end
 end
 match=match/length(subject.hdr.label);
 if match>best_match
  best_i=i;
  best_match=match;
 end
end

% now see if we have something
if best_match>0.9
 fprintf('Auto-detected this dataset as ''%s'', match=%2.1f%%\n',dataset_id(best_i).id,best_match*100)
 subject.detected_datatype=dataset_id(best_i).id;
 subject.dataset_mapping=all_m{best_i};
elseif isfield(subject,'repN') && subject.repN>10
 % likely an invasive dataset
 subject.dtype='EEG_invasive';
 if ~isfield(subject,'layout') || isempty(subject.layout)
  % try to auto-make a layout
  [subject.layout, subject.found_electrodes] = nmri_make_layout_invasive(subject);
  subject.valid_channels={'EEG'};
 end
 subject.montages=[]; % need to trust in eeg_score to deal with that
elseif best_match>0.7
 warning(sprintf('Detected this dataset as ''%s'', but only match=%2.1f%. Will try to continue with this%\n',dataset_id(best_i).id,best_match*100))
 subject.detected_datatype=dataset_id(best_i).id;
 subject.dataset_mapping=all_m{best_i};
else
 warning(sprintf('Could not reliably detect this dataset, best option is ''%s'', but only match=%2.1f%. Will continue with what we have in header.%\n',dataset_id(best_i).id,best_match*100))
 subject.detected_datatype='no match';
 % we may consider relying on the dataset - but hmmm
 if isfield(subject.hdr,'chantype')
  eeg_chan=sum(strcmpi(subject.hdr.chantype,'eeg'));
  meg_chan=sum(strcmpi(subject.hdr.chantype,'meg'));
  if meg_chan>0
   % if we have any MEG channels, probably this is an MEG dataset
   subject.dtype='MEG';
   % defaulting layout to CTF
   subject.layout='CTF275.lay';
   subject.montage='conf/montages/CTF275.mat';
   subject.elec_file='GSN-HydroCel-257.sfp';
   subject.elec_fiducials = {'FidNz','FidT9','FidT10'};
   
  elseif eeg_chan>10
   % otherwise probably EEG, at least 10 channels
   subject.dtype='EEG';
   subject.layout='EEG1010.lay'; %defaulting to 10-10
   subject.montage='conf/montages/Goe_Routine.mat';
  else
   error('Could not determine the dataset type, likely a problem with the fileformat or unusual channels.')
  end
 end
 subject.dataset_mapping=all_m{best_i};
end



% check our datatype first
if (~isfield(subject,'dtype') )
 % set from best match
 subject.dtype=dataset_id(best_i).dtype;
end

% layout
if (~isfield(subject,'layout') )
 % set from best match
 subject.layout=dataset_id(best_i).layout;
end

% montage
if (~isfield(subject,'montages') )
 % set from best match
 subject.montages=dataset_id(best_i).montage;
end


% elec
if (~isfield(subject,'elec_file') && isfield(dataset_id(best_i),'elec_file') )
 % set from best match
 subject.elec_file=dataset_id(best_i).elec_file;
end
if (~isfield(subject,'elec_fiducials') && isfield(dataset_id(best_i),'elec_fiducials') )
 % set from best match
 subject.elec_fiducials=dataset_id(best_i).elec_fiducials;
end


% now check channeltypes
for i=1:size(subject.hdr.label,1)
 match=find(strcmpi(subject.hdr.label(i,1),dataset_id(best_i).label));
 if ~isempty(match)
  if strcmpi(subject.hdr.chantype{i,1},'unknown')
  % try to determine from dataset_id if Fieldtrip failed 
   subject.hdr.chantype(i,1)=dataset_id(best_i).chantype(match);
   subject.hdr.chanunit(i,1)=dataset_id(best_i).chanunit(match);
  end
 
  % cover the case that a channel is not valid for some case, e.g. the Fpz
  % for Routine EEG in Göttingen
  if strcmpi(dataset_id(best_i).chantype(match),'unknown')
   % overwrite waht fieldtrip thinks
   subject.hdr.chantype(i,1)=dataset_id(best_i).chantype(match);
   subject.hdr.chanunit(i,1)=dataset_id(best_i).chanunit(match);
  end
  
 end
end




end

