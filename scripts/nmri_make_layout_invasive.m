function [ layout, found_elecs ] = nmri_make_layout_invasive( subject )
%[ layout ] = nmri_make_layout_invasive( subject )
%  Will try to auto-generate a layout for invasice SEEG data


%build a custom layout of all EEG channels
% needs some fieldtrip functions introduced with releases after 05/2019
if exist('ft_appendlayout','file')~=2
 error('Need a Fieldtrip version >= 20190527 to build an automatic invasive EEG layout. Stopping here.')
end

disp('Will try to generate a default layout for invasive EEG')

eeg=subject.hdr.label(strcmpi(subject.hdr.chantype,'eeg'));
found_elecs={};
for i=1:length(eeg)
 alpha=isstrprop(eeg{i},'alpha');
 if alpha(1)  %starts with a letter 
  % check if we know this one
  if ~any(strcmp(found_elecs,eeg{i}(alpha)))
   %find possible partners
   m=~cellfun(@isempty,regexp(eeg,['^' eeg{i}(alpha) '[0-9]*$']));
   if sum(m)>4 % smallest elec is 5 contacts for DIXI
    % check if we have all contacts - this woudl usually be the case for
    % invasive electrodes (unless some channels have been skipped
    this_labels=eeg(m);
    this_elec=eeg{i}(alpha);
    valid=1;
    for x=1:length(this_labels)
     if ~any(strcmp(this_labels,[this_elec num2str(x)]))
      % at least one is missing
      valid=0;
     end
    end
    if valid
     found_elecs{end+1}= this_elec;
    end
   end
  end
 end
end

% sort once alphabetically
found_elecs=sort(found_elecs);


% now build the electrode blocks
% make a symetric matrix for L/R channels
left=~cellfun(@isempty,regexp(found_elecs,['^L.*']));
right=~cellfun(@isempty,regexp(found_elecs,['^R.*']));
% now parse all left for right partners
bilat=false(1,length(found_elecs));
for i=1:length(left)
 if left(i)
  this=~cellfun(@isempty,regexp(found_elecs,['^R' found_elecs{i}(2:end) '[0-9]*$']));
  if any(this)
   bilat=(bilat|this);
   bilat(i)=true;
  end
 end
end

% now reorder the elecs with master logic
found_elecs=[found_elecs(bilat&left) found_elecs(bilat&right) found_elecs(~bilat)];

% and rebuidl  the logic
left=~cellfun(@isempty,regexp(found_elecs,['^L.*']));
right=~cellfun(@isempty,regexp(found_elecs,['^R.*']));
% now parse all left for right partners
bilat=false(1,length(found_elecs));
for i=1:length(left)
 if left(i)
  this=~cellfun(@isempty,regexp(found_elecs,['^R' found_elecs{i}(2:end) '[0-9]*$']));
  if any(this)
   bilat=(bilat|this);
   bilat(i)=true;
  end
 end
end


% now build a shaft for each electrode
thisL={};
for i=1:length(found_elecs)
 cfg=[];
 cfg.layout='vertical';
 cfg.direction='TB';
 cfg.width=0.2;
 cfg.height=0.05;
 cfg.channel=[found_elecs{i} '*'];
 thisL{end+1}=ft_prepare_layout(cfg,subject.hdr);
end

% now order the shafts
% bilateral to symetric - left first
layoutB=[];
layoutS=[];

if any(bilat)
 cfg=[];
 cfg.direction='horizontal';
 cfg.align='top';
 layoutL=ft_appendlayout(cfg,thisL{bilat&left});
 layoutR=ft_appendlayout(cfg,thisL{bilat&right});
 cfg=[];
 cfg.direction='horizontal';
 cfg.align='top';
 cfg.distance=0.1;
 layoutB=ft_appendlayout(cfg,layoutL,layoutR);
end

% now add the loners
if any(~bilat)
 layoutR=[];
 layoutL=[];
 cfg=[];
 cfg.direction='horizontal';
 cfg.align='top';
 if any((~bilat)&left)
  layoutL=ft_appendlayout(cfg,thisL{(~bilat)&left});
 end
 if any((~bilat)&right)
  layoutR=ft_appendlayout(cfg,thisL{(~bilat)&right});
 end
 
 if ~isempty(layoutL) && ~isempty(layoutR)
  cfg=[];
  cfg.direction='horizontal';
  cfg.align='top';
  cfg.distance=0.1;
  layoutS=ft_appendlayout(cfg,layoutL,layoutR);
 elseif ~isempty(layoutL)
  layoutS=layoutL;
 elseif ~isempty(layoutR)
  layoutS=layoutR;
 end
end

% now merge
if ~isempty(layoutS) && ~isempty(layoutB)
 cfg=[];
 cfg.direction='vertical';
 cfg.align='left';
 cfg.distance=0.1;
 layout=ft_appendlayout(cfg,layoutB,layoutS);
else
 if ~isempty(layoutS) 
  layout=layoutS;
 elseif ~isempty(layoutB) 
  layout=layoutB;
 else
  % probably no left-right oder
  cfg=[];
  cfg.direction='horizontal';
  cfg.align='top';
  layout=ft_appendlayout(cfg,thisL{:});
 end
end



end

