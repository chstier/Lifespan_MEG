function varargout = eeg_score(varargin)
% [data, subject]=eeg_score(cfg,data,subject)
%
% This script is both a EEG/MEG viewer and can be used to set events
% (e.g. "spikes") and do a trial-based rating for vigilance / sleep and
% techinal artifacts. It should be faster than ft_databrowser and somewhat
% more conventional vor EEG/MEG visual reading.
% Still, this is NOT a full replacement for a EEG/MEG clinical viewer!
%
% usage:
% cfg           = struct (array), can feature multiple montages, if wanted.
% 
%
% Last Modified by GUIDE v2.5 04-Oct-2019 16:26:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eeg_score_OpeningFcn, ...
                   'gui_OutputFcn',  @eeg_score_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before eeg_score is made visible.
function eeg_score_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eeg_score (see VARARGIN)

% Choose default command line output for eeg_score
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);





% check if we have cfg
if length(varargin)>0 && ~isempty(varargin{1})
 cfg=varargin{1};
else
 cfg=struct('empty',[]);
end
% check if we have data
if length(varargin)>1
 data=varargin{2};
else
 data=[];
end
% check if we have subject
if length(varargin)>2
 subject=varargin{3};
else
 subject=[];
end


% center main figure on screen
screensize = get(0, 'ScreenSize');
win_pos=get(hObject,'Position');
if isfield(cfg,'plot_layout') && cfg.plot_layout
 win_pos(1)=(screensize(3)-win_pos(3))/2; % center if layout
else
 win_pos(1)=(screensize(3)-win_pos(3))/3; % shifted a bit to the left
end
win_pos(2)=(screensize(4)-win_pos(4))/2;

set(hObject,'Position',win_pos);



if isfield(cfg,'vertical')
 if strcmpi(cfg.vertical,'spread')
  set(handles.vertical_align,'Value',1)
 elseif strcmpi(cfg.vertical,'overlay')
  set(handles.vertical_align,'Value',2)
 else
  % default to spread
  set(handles.vertical_align,'Value',1)
 end
end

if isfield(cfg,'hide_bad_channels')
 set(handles.hide_bad_channels,'Value',cfg.hide_bad_channels)
end

if isfield(cfg,'dash_bad_channels')
 set(handles.dash_bad_channels,'Value',cfg.dash_bad_channels)
end

% check data
if isempty(data)
 [filename, pathname] = uigetfile('*.mat');
 cfg.mfile=fullfile(pathname,filename);
end

if isfield(cfg,'mfile')
 if exist(cfg.mfile,'file')
  load(cfg.mfile)
  if  exist('data_source','var')
   % deal with source projected data
   data=data_source;
   clear data_source
   subject.hdr.label=data.label;
   subject.hdr.chanunit=repmat({'unknown'},size(data.label));
   subject.hdr.chantype=repmat({'unknown'},size(data.label));
  end
  if isempty(data) || ~isfield(data,'trial')
   error(['No data found in ' cfg.mfile])
  end
  data.mfile=cfg.mfile;
  % check for re-root
  if (~strcmp(cfg.mfile(1:length(subject.analysis_dir)),subject.analysis_dir))
   % there has been a move, now try to determine the true analysis dir
   analysis_dir=dirname(dirname(dirname(cfg.mfile))); %assmuming 3 level move
   subject=nmri_reroot_subject(subject,analysis_dir); 
  end
 else
  error(['File not found ' cfg.mfile])
 end
end

% check the dataset, if we have hdr info
if (~isfield(subject,'dtype') || ~isfield(subject,'montage')) && isfield(subject,'hdr')
 subject=nmri_determine_datatype(subject);
end

% set default cfg
if ~isfield(cfg,'montage') 
 if isfield(subject,'montage') && ~isempty(subject.montage) && exist(fullfile(subject.analysis_dir,subject.montage),'file')
  lcfg=load(fullfile(subject.analysis_dir,subject.montage),'cfg');
  montage=lcfg.cfg;
 elseif isfield(subject,'montages') && ~isempty(subject.montages) && exist(fullfile(subject.analysis_dir,subject.montages),'file')
  lcfg=load(fullfile(subject.analysis_dir,subject.montages),'cfg');
  montage=lcfg.cfg;
 elseif isfield(subject,'dtype') && strcmp(subject.dtype,'EEG_invasive') 
  % special behaviour for EEG invasive 
  montage.channel=data.label;
  montage.label=data.label;
  if isfield(subject,'found_electrodes')
   montage.color=repmat({'k'},size(montage.channel)); 
   alt={'b','k'};
   alt_i=1;
   for i=1:length(subject.found_electrodes)
    m=~cellfun(@isempty,regexp(montage.channel,['^' subject.found_electrodes{i} '[0-9]*$']));
    montage.color(m)=repmat(alt(alt_i),1,sum(m));
    alt_i=alt_i+1;
    if alt_i>length(alt)
     alt_i=1;
    end
   end
  else
   montage.color=repmat({'k'},size(cfg.channel));
  end
  montage.id='All Channels';
  montage.dtype='EEG';
  montage.def_scale=200;
  
  % make a 2:1 redux
  montage(2).id='Every 2nd channel';
  montage(2).label=montage(1).label(1:2:end);
  montage(2).channel=montage(1).channel(1:2:end);
  montage(2).color=montage(1).color(1:2:end);
  montage(2).dtype=montage(1).dtype;
  montage(2).def_scale=montage(1).def_scale;
  
  % make a 5:1 redux
  montage(3).id='Every 5th channel';
  montage(3).label=montage(1).label(1:5:end);
  montage(3).channel=montage(1).channel(1:5:end);
  montage(3).color=montage(1).color(1:5:end);
  montage(3).dtype=montage(1).dtype;
  montage(3).def_scale=montage(1).def_scale;
  
 else
  % not given, load all 
  if exist('subject') && isfield(subject,'analysis_dir')
   montage_dir=fullfile(subject.analysis_dir,'conf','montages');
  else
   % not set, default to dev dir
   montage_dir=fullfile('/data/projects/nfocke/dev','conf','montages');
  end
  if ~exist(montage_dir,'dir')
   error(['Could not find the montages in ' montage_dir])
  end
  all_d=dir(montage_dir);
  all_m={};
  for i=1:size(all_d,1)
   if ~strcmp(all_d(i).name(1),'.') && ~all_d(i).isdir
    all_m=[all_m {fullfile(montage_dir,all_d(i).name)}];
   end
  end
  montage=nf_load_mats_struct(all_m,'cfg');
 end
else
 % take from cfg
 montage=cfg.montage;
end


% now make bg_cols
if ~isfield(data,'trial_markings')
 
 % check for sampleinfo
 if ~isfield(data,'sampleinfo')
  w=warning('off','all'); 
  data=ft_checkdata(data, 'feedback', 'no', 'hassampleinfo', 'yes');
  warning(w);
 end
 
 data.trial_markings=cell(length(data.trial),4);
 % col1: sleep
 % col2: technical
 % col3: event
 % col4: rest/stimulation
 
 data.trial_markings_sampleinfo=cell(length(data.trial),2);
 % col1: sampleinfo start/stop
 % col2: seconds start/stop
 for i=1:length(data.trial)
  data.trial_markings_sampleinfo{i,1}=data.sampleinfo(i,:);
  data.trial_markings_sampleinfo{i,2}=data.sampleinfo(i,:)/data.fsample;
 end
end

if size(data.trial_markings,2)<4
 % make sure we have enough columns (backward compatabilty)
 for i=2:4
  if size(data.trial_markings,2)<i
   data.trial_markings=[data.trial_markings cell(size(data.trial_markings,1),1)];
  end
 end
end

bg_cols=make_bg_cols(data,subject);


% check for bad channels
if ~isfield(data,'bad_channels') 
 data.bad_channels={};
end

% register data
setappdata(handles.view_ax,'data',data);
setappdata(handles.view_ax,'cfg',cfg);
setappdata(handles.view_ax,'montage',montage);
setappdata(handles.view_ax,'active_trial',1);
setappdata(handles.view_ax,'trial_substep',1);
setappdata(handles.view_ax,'active_piece',[]);
setappdata(handles.view_ax,'xscale',1);
setappdata(handles.view_ax,'marker',[]);
setappdata(handles.view_ax,'marker_handle',[]);
setappdata(handles.view_ax,'event',[]);
setappdata(handles.view_ax,'event_handle',[]);
setappdata(handles.view_ax,'active_channel',[]); 

%set(handles.trial_slider,'Max',length(data.trial))
%set(handles.trial_slider,'Min',1)
%set(handles.trial_slider,'SliderStep',[1/(length(data.trial)-1) 1/(length(data.trial)-1)])
%set(handles.trial_slider,'Value',1)
set(handles.skip_next_technical,'Value',true);
set(handles.skip_next_vigilance,'Value',false);

% set montages
txt={};
for i=1:length(montage)
 if isfield(montage(i),'id')
  txt{end+1}=montage(i).id;
 else
  txt{end+1}=['Montage #' num2str(i)];
 end
end
set(handles.montage,'String',txt);
if isfield(cfg,'select_montage')
 if isnumeric(cfg.select_montage)
  set(handles.montage,'Value',cfg.select_montage);
 elseif ischar(cfg.select_montage)
  set(handles.montage,'Value',find(strcmpi(cfg.select_montage,txt)));
 end
end

% set yscale
if isfield(cfg,'def_scale')
 yscale=cfg.def_scale;
elseif isfield(montage(get(handles.montage,'Value')),'def_scale')
 yscale=montage(get(handles.montage,'Value')).def_scale;
else
 % not present, auto-determine
 yscale=range(reshape(data.trial{1},numel(data.trial{1}),1));
end
setappdata(handles.view_ax,'yscale',yscale);
set(handles.yscale,'String',num2str(yscale));

% set markings
if isfield(cfg,'score_vigilance')
 if cfg.score_vigilance==false
  % disable scoring of vigilance
  set(handles.vigilance,'Visible','Off','Enable','Off')
  set(handles.text_vigilance,'Visible','Off','Enable','Off')
  set(handles.skip_next_vigilance,'Visible','Off','Enable','Off')
 else
  % set auto-next for vigilance
  set(handles.skip_next_vigilance,'Value',true)
  set(handles.skip_next_technical,'Value',false)
 end
end
if isfield(cfg,'score_technical')
 if cfg.score_technical==false
  % disable scoring of vigilance
  set(handles.technical,'Visible','Off','Enable','Off')
  set(handles.text_technical,'Visible','Off','Enable','Off')
  set(handles.skip_next_technical,'Visible','Off','Enable','Off')
  % check if both are off, then can disable panel
  if strcmpi(get(handles.vigilance,'Visible'),'Off')
   set(handles.panel_ratings,'Visible','Off')
  end
 end
end
if isfield(cfg,'allow_events')
 if cfg.allow_events==false
  set(handles.panel_events,'Visible','Off')  
 end
else
 cfg.allow_events=true;
end
if ~isfield(cfg,'plot_layout') || ~isfield(subject,'layout')
 cfg.plot_layout=true;
end


% setup events/info listbox
%if cfg.allow_events
 event_list_show_Callback([],[],handles)
%end

% setup layout-channel
if cfg.plot_layout
 setappdata(handles.view_ax,'subject',subject)
 layout_show_Callback([],[],handles)
 subject=getappdata(handles.view_ax,'subject');
end


% set static subject info
txt={};
if isfield(data,'mfile')
 if length(data.mfile)>36;
  txt={sprintf('File: ...%s',data.mfile(end-36:end))};
 else
  txt={sprintf('File: %s',data.mfile)};
 end
else
 txt={'No File Info'};
end

if ~isempty(subject)
 txt{end+1}=sprintf('%s - %s (%s)',subject.id,subject.exam_id,subject.dtype);
else
 txt{end+1}='No Subject Info';
end

set(handles.subject_info,'String',txt);

% apply channel selection
select_channels(handles)

% make my colormap for bg and overview
mycolmap=[...
 1 1 1; ... % 1 - not rated
 0.8 1 0.8; ... % 2 - wake
 0.8 0.8 1; ... % 3 - S1
 0.7 0.7 1; ... % 4 - S2
 0.6 0.6 1; ... % 5 - S3
 0.8 0.6 1; ... % 6 - REM
 0.6 0.6 0.6; ... % 7 not used
 1 1 0.6; ... % 8 technical good, but not vigilance scored
 0.2 0.2 0; ... % 9 current trial
 1 0.5 0.8; ... % 10 - event
 1 0.5 0.5; ... % 11 - technical problem
 0.8 0.8 0.8; ... % 12 - eyes open
 1 0.7 0.7; ... % 13 - HV
 1 0.8 0.9; ... % 13 - PostHV
 1 1 0.5; ... % 14 - photostim
 
 ];
setappdata(handles.view_ax,'mycolmap',mycolmap);

setappdata(handles.view_ax,'bg_cols',bg_cols);
setappdata(handles.view_ax,'subject',subject)


%if cfg.allow_events
 update_evt_list(handles)
%end

% now draw the axis
update_view(handles)



% UIWAIT makes eeg_score wait for user response (see UIRESUME)
uiwait(handles.figure1);


function bg_cols=make_bg_cols(data,subject)
bg_cols=ones(1,length(data.trial));
% make bg's 
for i=1:length(data.trial)
 if (islogical(data.trial_markings{i,2}) &&~data.trial_markings{i,2})
  % technical problem - always red
  bg_cols(i)=11;
 elseif (ischar(data.trial_markings{i,4}) && ~strcmp(data.trial_markings{i,4},'eyes_closed') && ~strcmp(data.trial_markings{i,4},''))
  % not eyes closed, so give a color based on this
  if strcmp(data.trial_markings{i,4},'eyes_open') 
   bg_cols(i)=12;
  elseif strcmp(data.trial_markings{i,4},'HV')
   bg_cols(i)=13;
  elseif strcmp(data.trial_markings{i,4},'postHV')
   bg_cols(i)=14;
  elseif strcmp(data.trial_markings{i,4},'PS')
   bg_cols(i)=15;
  else
   bg_cols(i)=7;
  end  
 elseif strcmp(data.trial_markings{i,1},'wake') 
  % good
  bg_cols(i)=2;
 elseif isempty(data.trial_markings{i,1}) 
  %vigilance not rated
  if (islogical(data.trial_markings{i,2}) &&data.trial_markings{i,2})
   bg_cols(i)=8; % technically good, but vigilance not rated
  else
   bg_cols(i)=1;
  end
 elseif strcmp(data.trial_markings{i,1},'sleep1') 
  % s1
  bg_cols(i)=3;
 elseif strcmp(data.trial_markings{i,1},'sleep2') 
  % s2
  bg_cols(i)=4;
 elseif strcmp(data.trial_markings{i,1},'sleep3') 
  % s3 
  bg_cols(i)=5;
 elseif strcmp(data.trial_markings{i,1},'rem')  
  % REM 
  bg_cols(i)=6;
 else
  % we should actually never end here, but who knows
  bg_cols(i)=1;
 end
end
% check for events
if (isfield(subject,'evt_timings_seconds'))
 nSpikes=length(subject.evt_timings_seconds);
 for iTrial = 1:length(data.trial)
  for iSpike = 1:nSpikes
   if (data.time{iTrial}(1,end) > subject.evt_timings_seconds(iSpike)) && (data.time{iTrial}(1,1) < subject.evt_timings_seconds(iSpike))
    data.trial_markings{iTrial,3}='event';
    bg_cols(iTrial)=10;
   end
  end
 end
end



% --- Outputs from this function are returned to the command line.
function varargout = eeg_score_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=getappdata(handles.view_ax,'data');
varargout{2}=getappdata(handles.view_ax,'subject');
delete(hObject);



% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

data=getappdata(handles.view_ax,'data');
subject=getappdata(handles.view_ax,'subject');
cfg=getappdata(handles.view_ax,'cfg');
montage=getappdata(handles.view_ax,'montage');
active_trial=getappdata(handles.view_ax,'active_trial');
trial_substep=getappdata(handles.view_ax,'trial_substep');
active_piece=getappdata(handles.view_ax,'active_piece');
yscale=getappdata(handles.view_ax,'yscale');
xscale=getappdata(handles.view_ax,'xscale');
marker=getappdata(handles.view_ax,'marker');
event=getappdata(handles.view_ax,'event');
active_channel=getappdata(handles.view_ax,'active_channel');

do_update=false;

% check if a popup box has focus
if ~isempty(gco)
 try
  focus_obj=get(gco,'style');
 catch
  focus_obj='';
 end
else
 focus_obj='';
end

skip_next=false;
% master key switch
switch(eventdata.Key)
 case 'rightarrow'
  if ~isempty(eventdata.Modifier) && ~isempty(marker)
   % move marker
   if any(strcmp(eventdata.Modifier,'control'))
    marker=marker+(1/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'shift'))
    marker=marker+(5/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'alt'))
    marker=marker+(30/data.fsample);
   end
   if marker>data.time{active_trial}(active_piece(end))
    marker=data.time{active_trial}(active_piece(end));
   end
   setappdata(handles.view_ax,'marker',marker);
   update_marker(handles);
  elseif ~isempty(eventdata.Modifier) && ~isempty(event)
   % move event
   if any(strcmp(eventdata.Modifier,'control'))
    subject.evt_timings_seconds(event)=subject.evt_timings_seconds(event)+(1/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'shift'))
    subject.evt_timings_seconds(event)=subject.evt_timings_seconds(event)+(5/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'alt'))
    subject.evt_timings_seconds(event)=subject.evt_timings_seconds(event)+(30/data.fsample);
   end
   if subject.evt_timings_seconds(event)>data.time{active_trial}(active_piece(end))
    % go to next page
    trial_substep=trial_substep+1;
    if trial_substep>1/xscale
     active_trial=active_trial+1;
     trial_substep=1;
    end
    if active_trial>length(data.trial)
     active_trial=length(data.trial);
     trial_substep=1/xscale;
    end
    setappdata(handles.view_ax,'active_trial',active_trial);
    setappdata(handles.view_ax,'trial_substep',trial_substep);
    do_update=true;
   end
   % convert to samples, in some MEG datasets the time does not start with 0
   % find the right trial for this time, if any
   [m,i]=min(cellfun(@min,(cellfun(@(x) abs(x-subject.evt_timings_seconds(event)),data.time, 'UniformOutput', false))));
   if (m<0.5) % allow some margin for error  
    subject.evt_timings_sample(event)=data.sampleinfo(i,1)+round((subject.evt_timings_seconds(event)-data.time{i}(1))*data.fsample);
   else 
    % should not happen, but just in case
    subject.evt_timings_sample(event)=NaN;
   end

   setappdata(handles.view_ax,'subject',subject);
   update_event(handles); 
  else
   % go forward
   trial_substep=trial_substep+1;
   if trial_substep>1/xscale
    active_trial=active_trial+1;
    trial_substep=1;
   end
   if active_trial>length(data.trial)
    active_trial=length(data.trial);
    trial_substep=1/xscale;
   end
   setappdata(handles.view_ax,'active_trial',active_trial);
   setappdata(handles.view_ax,'trial_substep',trial_substep);
   setappdata(handles.view_ax,'event',[]);
   do_update=true;
  end
 case 'leftarrow'
  if ~isempty(eventdata.Modifier) && ~isempty(marker)
   % move marker
   if any(strcmp(eventdata.Modifier,'control'))
    marker=marker-(1/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'shift'))
    marker=marker-(5/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'alt'))
    marker=marker-(30/data.fsample);
   end
   if marker<data.time{active_trial}(active_piece(1))
    marker=data.time{active_trial}(active_piece(1));
   end
   setappdata(handles.view_ax,'marker',marker);
   update_marker(handles); 
  elseif ~isempty(eventdata.Modifier) && ~isempty(event)
   % move event
   if any(strcmp(eventdata.Modifier,'control')) 
    subject.evt_timings_seconds(event)=subject.evt_timings_seconds(event)-(1/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'shift'))
    subject.evt_timings_seconds(event)=subject.evt_timings_seconds(event)-(5/data.fsample);
   elseif all(strcmp(eventdata.Modifier,'alt'))
    subject.evt_timings_seconds(event)=subject.evt_timings_seconds(event)-(30/data.fsample);
   end
   if subject.evt_timings_seconds(event)<data.time{active_trial}(active_piece(1))
    % go back
    trial_substep=trial_substep-1;
    if trial_substep<=0
     active_trial=active_trial-1;
     trial_substep=1/xscale;
    end
    if active_trial<1
     active_trial=1;
     trial_substep=1;
    end
    setappdata(handles.view_ax,'active_trial',active_trial);
    setappdata(handles.view_ax,'trial_substep',trial_substep);
    do_update=true;
   end
   % convert to samples, in some MEG datasets the time does not start with 0
   % find the right trial for this time, if any
   [m,i]=min(cellfun(@min,(cellfun(@(x) abs(x-subject.evt_timings_seconds(event)),data.time, 'UniformOutput', false))));
   if (m<0.5) % allow some margin for error  
    subject.evt_timings_sample(event)=data.sampleinfo(i,1)+round((subject.evt_timings_seconds(event)-data.time{i}(1))*data.fsample);
   else 
    % should not happen, but just in case
    subject.evt_timings_sample(event)=NaN;
   end
   setappdata(handles.view_ax,'subject',subject);
   update_event(handles); 
  else
   % go back
   trial_substep=trial_substep-1;
   if trial_substep<=0
    active_trial=active_trial-1;
    trial_substep=1/xscale;
   end
   if active_trial<1
    active_trial=1;
    trial_substep=1;
   end
   setappdata(handles.view_ax,'active_trial',active_trial);
   setappdata(handles.view_ax,'trial_substep',trial_substep);
   setappdata(handles.view_ax,'event',[]);
   do_update=true;
  end

 case 'uparrow'
  if ~strcmp(focus_obj,'popupmenu')
   yscale=yscale*0.75;
   setappdata(handles.view_ax,'yscale',yscale);
   set(handles.yscale,'String',num2str(yscale));
   do_update=true;
  end
 case 'downarrow'
  if ~strcmp(focus_obj,'popupmenu')
   yscale=yscale*1.25;
   setappdata(handles.view_ax,'yscale',yscale);
   set(handles.yscale,'String',num2str(yscale));
   do_update=true;
  end
  
 case 'e'
 if ~isempty(marker) && strcmpi(get(handles.marking_set,'Enable'),'On')
  % call callback
  marking_set_Callback(handles.marking_set, [], handles)
 end
 
 case 'c'
  if ~isempty(active_channel)
   % call callback
   toggle_channel_Callback(handles.toggle_channel, [], handles)
  end
 
 
 case 'delete'
 if ~isempty(event) && strcmpi(get(handles.marking_remove,'Enable'),'On')
  % call callback
  marking_remove_Callback(handles.marking_remove, [], handles)
 end
 
end
% now deal with other special cases (multiple options)
if ~strcmp(focus_obj,'edit')
 
 % score vigilance
 if ~isfield(cfg,'score_vigilance') || cfg.score_vigilance
  if strcmp(eventdata.Key,'0') || strcmp(eventdata.Key,'numpad0')
   data.trial_markings{active_trial,1}='wake';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_vigilance,'Value');
  end
  if strcmp(eventdata.Key,'1') || strcmp(eventdata.Key,'numpad1')
   data.trial_markings{active_trial,1}='sleep1';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_vigilance,'Value');
  end
  if strcmp(eventdata.Key,'2') || strcmp(eventdata.Key,'numpad2')
   data.trial_markings{active_trial,1}='sleep2';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_vigilance,'Value');
  end
  if strcmp(eventdata.Key,'3') || strcmp(eventdata.Key,'numpad3')
   data.trial_markings{active_trial,1}='sleep3';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_vigilance,'Value');
  end
  if strcmp(eventdata.Key,'4') || strcmp(eventdata.Key,'numpad4')
   data.trial_markings{active_trial,1}='rem';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_vigilance,'Value');  
  end
 end
 
 % technical
 if ~isfield(cfg,'score_technical') || cfg.score_technical
  if strcmp(eventdata.Key,'return')
   data.trial_markings{active_trial,2}=true;
   setappdata(handles.view_ax,'data',data);
   skip_next=get(handles.skip_next_technical,'Value');
   do_update=true;
  end
  if strcmp(eventdata.Key,'backspace') || strcmp(eventdata.Key,'subtract') || strcmp(eventdata.Key,'hyphen')
   data.trial_markings{active_trial,2}=false;
   setappdata(handles.view_ax,'data',data);
   skip_next=get(handles.skip_next_technical,'Value');
   do_update=true;
  end
 end
 
 % stimulation / eyes / HV
 if ~isfield(cfg,'score_technical') || cfg.score_technical
  if strcmp(eventdata.Key,'7') || strcmp(eventdata.Key,'numpad7')
   data.trial_markings{active_trial,4}='eyes_closed';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_stimuli,'Value');
  end
  if strcmp(eventdata.Key,'8') || strcmp(eventdata.Key,'numpad8')
   data.trial_markings{active_trial,4}='eyes_open';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_stimuli,'Value');
  end
  if strcmp(eventdata.Key,'9') || strcmp(eventdata.Key,'numpad9')
   data.trial_markings{active_trial,4}='HV';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_stimuli,'Value');
  end
  if strcmp(eventdata.Key,'6') || strcmp(eventdata.Key,'numpad6')
   data.trial_markings{active_trial,4}='postHV';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_stimuli,'Value');
  end
  if strcmp(eventdata.Key,'comma')
   data.trial_markings{active_trial,4}='PS';
   setappdata(handles.view_ax,'data',data);
   do_update=true;
   skip_next=get(handles.skip_next_stimuli,'Value');
  end
 end
 
end

% now deal with skip
if skip_next
 update_view(handles)
 active_trial=active_trial+1;
 if active_trial>length(data.trial)
  active_trial=length(data.trial);
 end
 setappdata(handles.view_ax,'active_trial',active_trial);
 do_update=true;
end
if do_update
 % update bg_cols
 bg_cols=make_bg_cols(data,subject); 
 setappdata(handles.view_ax,'bg_cols',bg_cols);
 update_view(handles)
end

% --- Executes during object creation, after setting all properties.
function view_ax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate view_ax

function update_view(handles)
% does complete redraw
data=getappdata(handles.view_ax,'data');
cfg=getappdata(handles.view_ax,'cfg');
montage=getappdata(handles.view_ax,'montage');
active_trial=getappdata(handles.view_ax,'active_trial');
trial_substep=getappdata(handles.view_ax,'trial_substep');
yscale=getappdata(handles.view_ax,'yscale');
xscale=getappdata(handles.view_ax,'xscale');
chan_sel=getappdata(handles.view_ax,'chan_sel');
bg_cols=getappdata(handles.view_ax,'bg_cols');
mycolmap=getappdata(handles.view_ax,'mycolmap');
subject=getappdata(handles.view_ax,'subject');
marker=getappdata(handles.view_ax,'marker');
event=getappdata(handles.view_ax,'event');
active_channel=getappdata(handles.view_ax,'active_channel'); 

show_data=data.trial{active_trial}(chan_sel.sel,:);

% now start plotting
axes(handles.view_ax);
cla
setappdata(handles.view_ax,'marker_handle',[]);
hold on
if isempty(chan_sel.sel)
 % no channels
 text(0.35,-0.45,'no channels found for this montage')
 active_piece=[];
 active_time=[];
 spread=1;
else
 % plot data of active trial/substep
 active_piece=round(length(data.time{active_trial})*xscale*(trial_substep-1)+1):round(length(data.time{active_trial})*xscale*(trial_substep));
 active_time=data.time{active_trial}(active_piece);
 % check if we want to spread
 spread=get(handles.vertical_align,'Value');
 if spread==1
  % normalize with yscale and spread per channel
  show_data=(show_data/yscale);
  for i=1:length(chan_sel.sel)
   show_data(i,active_piece)=show_data(i,active_piece)+(length(chan_sel.sel)-i); % spread out channels
  end

  for i=1:length(chan_sel.sel)
   % check if active channel
   line_width=1;
   if ~isempty(active_channel) 
    if strcmpi(active_channel,data.label{chan_sel.sel(i)})
     line_width=2;
    end
   end
   % check for bad
   line_style='-';
   if isfield(data,'bad_channels') && get(handles.dash_bad_channels,'Value')==1 &&  any(strcmp(data.label{chan_sel.sel(i)},data.bad_channels))
    line_style=':';
   end
   % if have colors
   if isfield(chan_sel,'color') && ~isempty(chan_sel.color)
    plot(data.time{active_trial}(active_piece),show_data(i,active_piece),'Color',chan_sel.color{i},'HitTest','Off','LineWidth',line_width,'LineStyle',line_style); 
   else
    % no color 
    plot(data.time{active_trial}(active_piece),show_data(i,active_piece),'HitTest','Off','LineWidth',line_width,'LineStyle',line_style); 
   end
  end
  
 else
  % plot as is - no spread
  for i=1:length(chan_sel.sel)
   % check for bad
   line_style='-';
   if isfield(data,'bad_channels') && get(handles.dash_bad_channels,'Value')==1 && any(strcmp(data.label{chan_sel.sel(i)},data.bad_channels))
    line_style='-.';
   end
   if isfield(chan_sel,'color') && ~isempty(chan_sel.color) 
    plot(data.time{active_trial}(active_piece),show_data(i,active_piece),'Color',chan_sel.color{i},'HitTest','Off','LineStyle',line_style); 
   else
    % no color - plot as one
    plot(data.time{active_trial}(active_piece),show_data(i,active_piece),'HitTest','Off','LineStyle',line_style); 
   end
  end
 end
end
hold off

ax=gca;
if ~isempty(active_piece)
 xlim(ax,[data.time{active_trial}(active_piece(1)) data.time{active_trial}(active_piece(end))]);
end
if spread==1
 ylim(ax,[-1 length(chan_sel.sel)]); 
 % print labels if spread out
 chans=chan_sel.label([length(chan_sel.sel):-1:1]);
 %check for active channel
 if ~isempty(active_channel) 
  for i=1:length(chan_sel.sel)
   % check if active channel
   if strcmpi(active_channel,data.label{chan_sel.sel(i)})
    chans{length(chan_sel.sel)-i+1}=['\bf ' chans{length(chan_sel.sel)-i+1}];
   end
  end
 end
 
 %check for bad channel
 if isfield(data,'bad_channels') && ~isempty(data.bad_channels) 
  for i=1:length(chan_sel.sel)
   % check if bad
   if any(strcmp(data.bad_channels,data.label{chan_sel.sel(i)}))
    chans{length(chan_sel.sel)-i+1}=['\color{red} ' chans{length(chan_sel.sel)-i+1} ];
   end
  end
 end
 
 set(ax,'YTick',[1:length(chan_sel.sel)]-1,'YTickLabel',chans)
else
 ylim(ax,[-yscale yscale])
 set(ax,'YTickMode','auto','YTickLabelMode','auto') 
end
set(ax,'XMinorGrid','On')
%set(handles.trial_slider,'Value',active_trial)

%now update markings
if isempty(data.trial_markings{active_trial,1})
 sel=6;
else
 switch data.trial_markings{active_trial,1}
  case 'wake'
   sel=1;
  case 'sleep1'
   sel=2;
  case 'sleep2'
   sel=3;
  case 'sleep3'
   sel=4;
  case 'rem'
   sel=5;
  otherwise
   sel=6;
 end
end
set(handles.vigilance,'Value',sel);
if isempty(data.trial_markings{active_trial,2})
 tsel=3;
 set(ax,'Color',[1 1 1])
elseif data.trial_markings{active_trial,2}
 tsel=1;
 set(ax,'Color',[1 1 1])
else
 tsel=2;
end
set(handles.technical,'Value',tsel);

if isempty(data.trial_markings{active_trial,4})
 ssel=6;
elseif strcmp(data.trial_markings{active_trial,4},'eyes_closed')
 ssel=1;
elseif strcmp(data.trial_markings{active_trial,4},'eyes_open')
 ssel=2;
elseif strcmp(data.trial_markings{active_trial,4},'HV')
 ssel=3;
elseif strcmp(data.trial_markings{active_trial,4},'postHV')
 ssel=4;
elseif strcmp(data.trial_markings{active_trial,4},'PS')
 ssel=5;
else
 ssel=6;
end
set(handles.stimuli,'Value',ssel);



% check for event
% event stays alltime
% check which one it is
if isfield(subject,'evt_timings_seconds')
 pyl=ylim(ax); 
 nSpikes=length(subject.evt_timings_seconds);

 %check delete button
 evt_visible=0;
 
 % parse events
 for iSpike = 1:nSpikes
  if (active_time(1,end) > subject.evt_timings_seconds(iSpike)) && (active_time(1,1) < subject.evt_timings_seconds(iSpike))
  
   % plot and register callback
   this_evt=line(ax,[subject.evt_timings_seconds(iSpike) subject.evt_timings_seconds(iSpike)],[pyl(1) pyl(2)],'Color','r','LineWidth',0.5,'LineStyle','-.',...
   'ButtonDownFcn',@(hObject,eventdata)eeg_score('event_selection_Callback',hObject,eventdata,guidata(hObject)),'UserData',iSpike);
   % check if active event and save handle
   if event==iSpike
    setappdata(handles.view_ax,'event_handle',this_evt);
    set(this_evt,'Color',[0.8 0 0.8],'LineWidth',1.5)
    update_event_controls(handles)
    evt_visible=1;
   end
   % add ID is set
   if isfield(subject,'evt_IDs') && ischar(subject.evt_IDs{iSpike})
    text(ax,subject.evt_timings_seconds(iSpike),pyl(2),subject.evt_IDs{iSpike},'Tag',['Evt_ID_' num2str(iSpike)])
   end
  end
 end
end

if exist('evt_visible','var') && evt_visible==1
 set(handles.marking_remove,'Enable','On')
else
 set(handles.marking_remove,'Enable','Off')
end

if isfield(subject,'info_timings_seconds')
 nInfos=length(subject.info_timings_seconds);

 % parse infos
 for iInfos = 1:nInfos
  if (active_time(1,end) > subject.info_timings_seconds(iInfos)) && (active_time(1,1) < subject.info_timings_seconds(iInfos))
   % plot and register callback
   line(ax,[subject.info_timings_seconds(iInfos) subject.info_timings_seconds(iInfos)],[pyl(1) pyl(2)],'Color','m','LineWidth',0.5,'LineStyle',':');
   % add ID is set
   if isfield(subject,'info_IDs') && ischar(subject.info_IDs{iInfos})
    text(ax,subject.info_timings_seconds(iInfos),pyl(2),subject.info_IDs{iInfos},'Tag',['Info_ID_' num2str(iInfos)])
   end
  end
 end
end

% if (tsel==2)
%  % technical problem - always red
%  bg_cols(active_trial)=11;
% elseif sel==1 
%  % good
%  bg_cols(active_trial)=2;
% elseif sel==6
%  %vigilance not rated
%  if tsel==1
%   bg_cols(active_trial)=8;
%  else
%   bg_cols(active_trial)=1;
%  end
% elseif sel==2
%  % s1
%  bg_cols(active_trial)=3;
% elseif sel==3
%  % s2
%  bg_cols(active_trial)=4;
% elseif sel==4
%  % s3 
%  bg_cols(active_trial)=5;
% elseif sel==5 
%  % REM 
%  bg_cols(active_trial)=6;
% end
if bg_cols(active_trial)~=10
 set(ax,'Color',mycolmap(bg_cols(active_trial),:))
else
 % for event we do differently
 if (tsel==2)
  this_col=11;
 elseif sel==1 
  this_col=2;
 elseif sel==6
  this_col=1;
 elseif sel==2
  this_col=3;
 elseif sel==3
  this_col=4;
 elseif sel==4
  this_col=5;
 else 
  this_col=6;
 end
 set(ax,'Color',mycolmap(this_col,:))
end

% now draw marker if set
if ~isempty(marker)
 % check if on current page
 if ~isempty(active_piece)
  if marker>=data.time{active_trial}(active_piece(1)) && marker<=data.time{active_trial}(active_piece(end))
   update_marker(handles);
   % enable set at marker
   set(handles.marking_set,'Enable','On')
   set(handles.marking_remove,'Enable','Off')
   update_event_controls(handles)
  else
   set(handles.marking_set,'Enable','Off')
  end
 end
end

% save active piece
setappdata(handles.view_ax,'active_piece',active_piece);

% now update overview
setappdata(handles.view_ax,'bg_cols',bg_cols);
% mark current trial
bg_cols(active_trial)=9;
im=image(handles.overview,bg_cols);
colormap(handles.overview,mycolmap)
set(handles.overview,'FontSize',8,'YTick',[],'YTickLabel',[])
set(im,'ButtonDownFcn',@(hObject,eventdata)eeg_score('overview_ButtonDownFcn',hObject,eventdata,guidata(hObject)))

% update info text - trials
txt={sprintf('Trial: %d/%d (%d), ToDo: %d',active_trial,length(data.trial),trial_substep,sum(cellfun(@isempty,data.trial_markings(:,1)))),...
sprintf('Wake: %d, Sleep: %d, Bad: %d',sum(strcmp(data.trial_markings(:,1),'wake')),(sum(strcmp(data.trial_markings(:,1),'sleep1'))...
+sum(strcmp(data.trial_markings(:,1),'sleep2'))+sum(strcmp(data.trial_markings(:,1),'sleep3'))+sum(strcmp(data.trial_markings(:,1),'rem'))),...
sum(cellfun(@(x) islogical(x) && x==false,data.trial_markings(:,2))))};
set(handles.info_txt,'String',txt)

% update info text - events
if isfield(subject,'evt_timings_seconds')
 evt_c=num2str(length(subject.evt_timings_seconds));
 if ~isempty(event)
  evt_act=sprintf('#%d, ID=%s @ %5.4f',event,subject.evt_IDs{event},subject.evt_timings_seconds(event));
 else
  evt_act='none';
 end
else
 evt_c='none';
 evt_act='none';
end
txt={sprintf('Events: %s, ActiveEvent: %s',evt_c,evt_act)};
set(handles.event_info,'String',txt)

% update channel selection
txt={};
if isempty(active_channel) 
 txt={'No channel selected'};
else
 idx=find(strcmp(active_channel,data.label));
 if ~isempty(idx)
  label=chan_sel.label(find(chan_sel.sel==idx,1));
  if ~isempty(label)
   txt={sprintf('Label: %s, Channel: %s',label{1},active_channel)};
  else
   txt={sprintf('Channel: %s',active_channel)};
  end
  status='Good';
  if isfield(data,'bad_channels') && any(strcmp(active_channel,data.bad_channels))
   status='Bad';
  end
  txt=[txt { ['Status: ' status] }]; 
 else
  txt={sprintf('Channel: %s, NOT IN DATA',active_channel)};
 end
end
set(handles.channel_info,'String',txt);

% check event controls
update_event_controls(handles)



function update_evt_list(handles)
data=getappdata(handles.view_ax,'data');
subject=getappdata(handles.view_ax,'subject');
marker=getappdata(handles.view_ax,'marker');
event=getappdata(handles.view_ax,'event');

% udpate event window, if any
if isappdata(handles.view_ax,'evFig')
 evFig=getappdata(handles.view_ax,'evFig');
 evLb=getappdata(handles.view_ax,'evLb');
 if ishandle(evLb)
  % make a new list of all events and info times and sort unique
  conc_evts_time=[];
  % cont here, check if all is set
  if isfield(subject,'evt_timings_seconds')
   if isfield(subject,'info_timings_seconds')
    conc_evts_time=unique([subject.evt_timings_seconds;subject.info_timings_seconds]);
   else
    conc_evts_time=unique(subject.evt_timings_seconds);
   end
  elseif isfield(subject,'info_timings_seconds')
   conc_evts_time=unique(subject.info_timings_seconds);
  end
  
  % now go through and build txt and jumplist
  txt={};
  jumplist=[];
  for i=1:length(conc_evts_time)
   if isfield(subject,'evt_timings_seconds')
    evtf=find(subject.evt_timings_seconds==conc_evts_time(i));
    for ii=1:length(evtf) 
     secs=subject.evt_timings_seconds(evtf(ii));
     % check if there is a trial with this event in the data
     found_trial=[];
     if secs>=0 && secs<=data.time{end}(1,end)
      % need to parse
      for ti=1:length(data.trial)
       if secs>=data.time{ti}(1,1) && secs<=data.time{ti}(1,end)
        found_trial=ti;
        continue
       end
      end
     end
     h=floor(secs/3600);
     secs=secs-(h*3600);
     m=floor(secs/60);
     secs=secs-(m*60);    
     if ~isempty(found_trial)     
      txt=[txt;[sprintf('%02.0f:%02.0f:%02.0f %s [EVT]',h,m,secs,subject.evt_IDs{evtf(ii)})] ];
      jumplist=[jumplist ; [conc_evts_time(i) evtf(ii) found_trial]];
     else
      txt=[txt;[sprintf('(%02.0f:%02.0f:%02.0f %s [EVT])',h,m,secs,subject.evt_IDs{evtf(ii)})] ];
      jumplist=[jumplist ; [conc_evts_time(i) NaN NaN]]; 
     end
    end
   end
  
   if isfield(subject,'info_timings_seconds')
    infof=find(subject.info_timings_seconds==conc_evts_time(i));
    for ii=1:length(infof) 
     secs=subject.info_timings_seconds(infof(ii));
     % check if there is a trial with this event in the data
     found_trial=[];
     if secs>=0 && secs<=data.time{end}(1,end)
      % need to parse
      for ti=1:length(data.trial)
       if secs>=data.time{ti}(1,1) && secs<=data.time{ti}(1,end)
        found_trial=ti;
        continue
       end
      end
     end
     h=floor(secs/3600);
     secs=secs-(h*3600);
     m=floor(secs/60);
     secs=secs-(m*60);
    
     if ~isempty(found_trial) 
      txt=[txt;[sprintf('%02.0f:%02.0f:%02.0f %s [INFO]',h,m,secs,subject.info_IDs{infof(ii)})] ];
      jumplist=[jumplist ; [conc_evts_time(i) NaN found_trial]];
     else
      txt=[txt;[sprintf('(%02.0f:%02.0f:%02.0f %s [INFO])',h,m,secs,subject.info_IDs{infof(ii)})] ];
      jumplist=[jumplist ; [conc_evts_time(i) NaN NaN]];
     end
    end
   end
  end
  if isempty(txt)
   txt={'No markers'};
  end

  set(evLb,'String',txt)
  setappdata(handles.view_ax,'evJumplist',jumplist);
 end
end


function update_marker(handles)
marker=getappdata(handles.view_ax,'marker');
marker_handle=getappdata(handles.view_ax,'marker_handle');
% now draw marker if set
if ~isempty(marker)
 pyl=ylim(handles.view_ax); 
 % check if set
 if isempty(marker_handle)
  % draw
  marker_handle=line(handles.view_ax,[marker marker],[pyl(1) pyl(2)],'Color','b','LineWidth',0.5,'LineStyle','--');
  setappdata(handles.view_ax,'marker_handle',marker_handle);
 else
  % move
  set(marker_handle,'XData',[marker marker])
 end
end

function update_event(handles)
event=getappdata(handles.view_ax,'event');
event_handle=getappdata(handles.view_ax,'event_handle');
subject=getappdata(handles.view_ax,'subject');
data=getappdata(handles.view_ax,'data');
active_trial=getappdata(handles.view_ax,'active_trial');
% now draw marker if set
if ~isempty(event)
 % check if set
 if ~isempty(event_handle)
  spikes=round(subject.evt_timings_seconds*data.fsample);
  ti=spikes(event)-data.sampleinfo(active_trial,1);
  % make active
  try
   set(event_handle,'XData',[data.time{active_trial}(ti) data.time{active_trial}(ti)],'Color',[0.8 0 0.8], 'LineWidth', 1.5)
   % check text
   if isfield(subject,'evt_IDs')
    txt_handle=findobj('Tag',['Evt_ID_' num2str(event)]);
    if ~isempty(txt_handle)
     pos=get(txt_handle,'Position');
     set(txt_handle,'Position',[data.time{active_trial}(ti)  pos(2)])
    end
   end
  catch
   % seems to be deleted
   setappdata(handles.view_ax,'event_handle',[]);
   setappdata(handles.view_ax,'event',[]);
   event=[];
   set(handles.marking_remove,'Enable','Off')
   set(handles.event_name,'Enable','Off')
   set(handles.event_name_custom,'Enable','Off')
   update_event_controls(handles)
  end
 end
 % update info text - events
 if isfield(subject,'evt_timings_seconds')
  evt_c=num2str(length(subject.evt_timings_seconds));
  if ~isempty(event)
   evt_act=sprintf('#%d, ID=%s @ %5.4f',event,subject.evt_IDs{event},subject.evt_timings_seconds(event));
  else
   evt_act='none';
  end
 else
  evt_c='none';
  evt_act='none';
 end
 txt={sprintf('Events: %s, ActiveEvent: %s',evt_c,evt_act)};
 set(handles.event_info,'String',txt)
 
  %deal with event window
  if isappdata(handles.view_ax,'evFig')
   evFig=getappdata(handles.view_ax,'evFig');
   evLb=getappdata(handles.view_ax,'evLb');
   if ishandle(evLb)
    % still active = not closed, so select active
    jumplist=getappdata(handles.view_ax,'evJumplist');
    fi=find(jumplist(:,2)==event);
    if ~isempty(fi)
     set(evLb,'Value',fi);
    end
   end
  end
 
end

function update_event_controls(handles)
% will update event controls based in selection
event=getappdata(handles.view_ax,'event');
subject=getappdata(handles.view_ax,'subject');
if ~isempty(event)
 if subject.evt_timings_seconds(event)<max(subject.evt_timings_seconds)
  set(handles.event_forward,'Enable','On')
  set(handles.event_last,'Enable','On')
 else
  set(handles.event_forward,'Enable','Off')
  set(handles.event_last,'Enable','Off')
 end
 if subject.evt_timings_seconds(event)>min(subject.evt_timings_seconds)
  set(handles.event_backward,'Enable','On')
  set(handles.event_first,'Enable','On')
 else
  set(handles.event_backward,'Enable','Off')
  set(handles.event_first,'Enable','Off')
 end
else
 % no event selected
 if isfield(subject,'evt_timings_seconds')
  set(handles.event_first,'Enable','On')
  set(handles.event_last,'Enable','On')
 else
  set(handles.event_first,'Enable','Off')
  set(handles.event_last,'Enable','Off')
 end
 set(handles.event_backward,'Enable','Off')
 set(handles.event_forward,'Enable','Off')
end


function select_channels(handles)
% this will check for excluded channels
data=getappdata(handles.view_ax,'data');
cfg=getappdata(handles.view_ax,'cfg');
montage=getappdata(handles.view_ax,'montage');
sel_m=get(handles.montage,'value');
count_c=1;
for i=1:length(montage(sel_m).channel)
 
 % skip bad channels, if so requested and data
 if get(handles.hide_bad_channels,'Value')==1 && isfield(data,'bad_channels') && ~isempty(data.bad_channels) && any(strcmp(montage(sel_m).channel(i),data.bad_channels))
  continue
 end
 
 f=find(strcmpi(data.label,montage(sel_m).channel{i}));
 if ~isempty(f)
  chan_sel.sel(count_c,1)=f;
  chan_sel.label{count_c,1}=montage(sel_m).label{i};
  if ~isempty(montage(sel_m).color)
   chan_sel.color{count_c,1}=montage(sel_m).color{i};
  end
  count_c=count_c+1;
 end
end
if ~exist('chan_sel','var')
 % seems no channels found
 chan_sel.sel=[];
 chan_sel.label={};
 chan_sel.color={};
else
 % now set the right space for main plot
 maxl=max(cellfun('length',chan_sel.label));
 % read mainplot pos
 old_pos=get(handles.view_ax,'Position');
 % set new x pos based on max
 win_pos=get(handles.figure1,'Position');
 shiftx=(maxl*6)+20;
 set(handles.view_ax,'Position',[shiftx,old_pos(2),win_pos(3)-shiftx-10,old_pos(4)]);
end
setappdata(handles.view_ax,'chan_sel',chan_sel);

% --- Executes on selection change in vigilance.
function vigilance_Callback(hObject, eventdata, handles)
% hObject    handle to vigilance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vigilance contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vigilance
data=getappdata(handles.view_ax,'data');
subject=getappdata(handles.view_ax,'subject');
active_trial=getappdata(handles.view_ax,'active_trial');
switch get(hObject,'Value')
 case 1
  data.trial_markings{active_trial,1}='wake';
 case 2
  data.trial_markings{active_trial,1}='sleep1';
 case 3
  data.trial_markings{active_trial,1}='sleep2'; 
 case 4
  data.trial_markings{active_trial,1}='sleep3';
 case 5
  data.trial_markings{active_trial,1}='rem';
 case 6
  data.trial_markings{active_trial,1}=[];
end
setappdata(handles.view_ax,'data',data);
% update bg_cols
bg_cols=make_bg_cols(data,subject); 
setappdata(handles.view_ax,'bg_cols',bg_cols);
update_view(handles)


% --- Executes during object creation, after setting all properties.
function vigilance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vigilance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in technical.
function technical_Callback(hObject, eventdata, handles)
% hObject    handle to technical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns technical contents as cell array
%        contents{get(hObject,'Value')} returns selected item from technical
data=getappdata(handles.view_ax,'data');
subject=getappdata(handles.view_ax,'subject');
active_trial=getappdata(handles.view_ax,'active_trial');
switch get(hObject,'Value')
 case 1
  data.trial_markings{active_trial,2}=true;
 case 2
  data.trial_markings{active_trial,2}=false;
 case 3
  data.trial_markings{active_trial,2}=[];
end
setappdata(handles.view_ax,'data',data);
% update bg_cols
bg_cols=make_bg_cols(data,subject); 
setappdata(handles.view_ax,'bg_cols',bg_cols);
update_view(handles)

% --- Executes during object creation, after setting all properties.
function technical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to technical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
 % The GUI is still in UIWAIT, us UIRESUME
 uiresume(hObject);
else
 % The GUI is no longer waiting, just close it
 delete(hObject);
end
%deal with event window
if isappdata(handles.view_ax,'evFig')
 evFig=getappdata(handles.view_ax,'evFig');
 if ishandle(evFig)
  close(evFig)
 end
end
%deal with layout window
if isappdata(handles.view_ax,'layFig')
 layFig=getappdata(handles.view_ax,'layFig');
 if ishandle(layFig)
  close(layFig)
 end
end




% --- Executes on slider movement.
function trial_slider_Callback(hObject, eventdata, handles)
% hObject    handle to trial_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
active_trial=round(get(hObject,'Value'));
setappdata(handles.view_ax,'active_trial',active_trial);
update_view(handles)


% --- Executes during object creation, after setting all properties.
function trial_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in skip_next_technical.
function skip_next_technical_Callback(hObject, eventdata, handles)
% hObject    handle to skip_next_technical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skip_next_technical


% --- Executes on selection change in montage.
function montage_Callback(hObject, eventdata, handles)
% hObject    handle to montage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns montage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from montage

% update channel selection
select_channels(handles)
% redraw
update_view(handles)


% --- Executes during object creation, after setting all properties.
function montage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to montage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in marking_set.
function marking_set_Callback(hObject, eventdata, handles)
% hObject    handle to marking_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set a new event at marker postion, if not already present
event=getappdata(handles.view_ax,'event');
marker=getappdata(handles.view_ax,'marker');
subject=getappdata(handles.view_ax,'subject');
data=getappdata(handles.view_ax,'data');
bg_cols=getappdata(handles.view_ax,'bg_cols');

if ~isempty(marker)
 if ~isfield(subject,'evt_timings_seconds')
  subject.evt_timings_seconds=[];
  subject.evt_timings_sample=[];
  subject.evt_IDs={};
 end
 if ~isfield(subject,'evt_IDs')
  subject.evt_IDs=cell(length(subject.evt_timings_seconds),1);
 end
 n_evt=length(subject.evt_timings_seconds);
 subject.evt_timings_seconds(n_evt+1,1)=marker;
 % convert to samples, in some MEG datasets the time does not start with 0
 % find the right trial for this time, if any
 [m,i]=min(cellfun(@min,(cellfun(@(x) abs(x-marker),data.time, 'UniformOutput', false))));
 if (m<0.5) % allow some margin for error  
  subject.evt_timings_sample(n_evt+1,1)=data.sampleinfo(i,1)+round((marker-data.time{i}(1))*data.fsample);
 else 
  % should not happen, but just in case
  subject.evt_timings_sample(n_evt+1,1)=NaN;
 end
 
 contents = cellstr(get(handles.event_name,'String'));
 this_id=contents{get(handles.event_name,'Value')};
 if (strcmp(this_id,'...custom'))
  this_id=get(handles.event_name_custom,'String');
 end
 subject.evt_IDs{n_evt+1,1}=this_id; 
 setappdata(handles.view_ax,'subject',subject);
 % now unset marker and focus event
 setappdata(handles.view_ax,'marker',[]);
 setappdata(handles.view_ax,'event',n_evt+1);
 % update bg_cols
 bg_cols=make_bg_cols(data,subject); 
 setappdata(handles.view_ax,'bg_cols',bg_cols);
 % set buttons
 set(handles.marking_remove,'Enable','On')
 set(handles.marking_set,'Enable','Off')
 update_evt_list(handles)
 update_event(handles)
 update_event_controls(handles)

 % redraw
 update_view(handles)
end


% --- Executes on button press in marking_remove.
function marking_remove_Callback(hObject, eventdata, handles)
% hObject    handle to marking_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
event=getappdata(handles.view_ax,'event');
marker=getappdata(handles.view_ax,'marker');
subject=getappdata(handles.view_ax,'subject');
data=getappdata(handles.view_ax,'data');
if ~isempty(event) 
 % user confirm
 button=questdlg({'Do you really want to remove the active event?',['ID='   subject.evt_IDs{event} ', Time=' num2str(subject.evt_timings_seconds(event))]},'Confirm Event Deletion','No');
 if strcmpi(button,'Yes')
  % do it
  subject.evt_timings_seconds(event)=[];
  subject.evt_timings_sample(event)=[];
  subject.evt_IDs(event)=[];
  setappdata(handles.view_ax,'subject',subject);
  setappdata(handles.view_ax,'event',[]);
  % update bg_cols
  bg_cols=make_bg_cols(data,subject); 
  setappdata(handles.view_ax,'bg_cols',bg_cols);
  update_view(handles)
 end
end

function yscale_Callback(hObject, eventdata, handles)
% hObject    handle to yscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yscale as text
%        str2double(get(hObject,'String')) returns contents of yscale as a double
setappdata(handles.view_ax,'yscale',str2double(get(hObject,'String')));
% redraw
update_view(handles)


% --- Executes during object creation, after setting all properties.
function yscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xscale.
function xscale_Callback(hObject, eventdata, handles)
% hObject    handle to xscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xscale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xscale
contents = cellstr(get(hObject,'String'));
xscale=str2double(contents{get(hObject,'Value')});
trial_substep=getappdata(handles.view_ax,'trial_substep');
event=getappdata(handles.view_ax,'event');
old_xscale=getappdata(handles.view_ax,'xscale');
subject=getappdata(handles.view_ax,'subject');
data=getappdata(handles.view_ax,'data');
active_trial=getappdata(handles.view_ax,'active_trial');
active_piece=getappdata(handles.view_ax,'active_piece');
marker=getappdata(handles.view_ax,'marker');
setappdata(handles.view_ax,'xscale',xscale);
% check if we want to keep focus (marker or event)
focus=[];
if ~isempty(event)
 focus=subject.evt_timings_seconds(event);
elseif ~isempty(marker)
 focus=marker;
else
 % keep left aligned
 focus=data.time{active_trial}(active_piece(1));
end
if ~isempty(focus)
 % keep this in sight
 trial_substep=floor(nearest(data.time{active_trial},focus)/(length(data.time{active_trial})*xscale))+1;
end
if 1/trial_substep<xscale
 % check if we get too high
 trial_substep=1/xscale;
end
setappdata(handles.view_ax,'trial_substep',trial_substep);
% redraw
update_view(handles)


% --- Executes during object creation, after setting all properties.
function xscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function view_ax_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to view_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subject=getappdata(handles.view_ax,'subject');
cfg=getappdata(handles.view_ax,'cfg');
chan_sel=getappdata(handles.view_ax,'chan_sel');
data=getappdata(handles.view_ax,'data');
active_trial=getappdata(handles.view_ax,'active_trial');
active_piece=getappdata(handles.view_ax,'active_piece');
yscale=getappdata(handles.view_ax,'yscale');
active_channel=getappdata(handles.view_ax,'active_channel'); 

% left mouse click
if strcmp(eventdata.EventName,'Hit') && eventdata.Button==1 
 
 % check, if label or plot was clicked
 if ~isempty(active_trial) && ~isempty(active_piece) 
  ref=data.time{active_trial}(active_piece(1));
  if (eventdata.IntersectionPoint(1)<ref)
   % clicked on label
   spread=get(handles.vertical_align,'Value');
   if spread==1
    % and spread out
    chan=length(chan_sel.sel)-(round(eventdata.IntersectionPoint(2)));
    if chan>0 && chan<=length(chan_sel.sel)
     fprintf('Clicked on label=%s, channel=%s\n',chan_sel.label{chan},data.label{chan_sel.sel(chan)})
     if isempty(active_channel) || ~strcmp(active_channel,data.label{chan_sel.sel(chan)})
      setappdata(handles.view_ax,'active_channel',data.label{chan_sel.sel(chan)}); 
      set(handles.toggle_channel,'Enable','On')
     else
      setappdata(handles.view_ax,'active_channel',[]);
      set(handles.toggle_channel,'Enable','Off')
     end
     % call full redraw
     update_view(handles)
     % update layout
     if isappdata(handles.view_ax,'layFig')
      layFig=getappdata(handles.view_ax,'layFig');
      if ishandle(layFig)
       update_layout_channels(handles)
      end
     end
    end
   end  
   
  elseif (~isfield(cfg,'allow_events') || cfg.allow_events)
   % clicked on plot, reposition line marker if we allow events

   setappdata(handles.view_ax,'marker',eventdata.IntersectionPoint(1)); % set to X - time
   % uncheck event if set
   event_handle=getappdata(handles.view_ax,'event_handle');
   if ~isempty(event_handle)
    % set to normal again
    try
     set(event_handle,'Color','r','LineWidth',0.5)
    catch
     % may have been deleted - we do not care
    end
   end
   % unset event
   setappdata(handles.view_ax,'event_handle',[]);
   setappdata(handles.view_ax,'event',[]);
   % set contrls
   set(handles.marking_set,'Enable','On')
   set(handles.marking_remove,'Enable','Off')
   update_event_controls(handles)
   % popuplate name suggestions - retaining the active one, just in case
   contents = cellstr(get(handles.event_name,'String'));
   old_val=contents{get(handles.event_name,'Value')};
   % build new
   evt_IDs={'Spike 1';'Spike 2';'Spike 3'};
   if isfield(subject,'evt_IDs')
    evt_IDs=sort(unique(vertcat(evt_IDs,subject.evt_IDs)));
   end
   % add fixed at the beginning or end
   if ~any(strcmpi('',evt_IDs))
    evt_IDs=vertcat({''},evt_IDs);
   end
   if ~any(strcmpi('Unclear',evt_IDs))
    evt_IDs=vertcat(evt_IDs,{'Unclear'});
   end
   if ~any(strcmpi('...custom',evt_IDs))
    evt_IDs=vertcat(evt_IDs,{'...custom'});
   end
   % add IDs
   set(handles.event_name,'String',evt_IDs)
   if ~isempty(old_val)
    selected=find(strcmpi(old_val,evt_IDs));
    if isempty(selected)
     selected=2;
    end
    set(handles.event_name,'Value',selected)
   else
    set(handles.event_name,'Value',2)
   end
   set(handles.event_name,'Enable','On')
   % check custom
   if get(handles.event_name,'Value')==length(evt_IDs)
    set(handles.event_name_custom,'Enable','On')
   else
    set(handles.event_name_custom,'Enable','Off') 
   end
   % call small redraw
   update_marker(handles)
  end
 end
end

function event_selection_Callback(hObject,eventdata,handles)
% Callback for click on event marker - do highlight
% check for other 
event=getappdata(handles.view_ax,'event');
event_handle=getappdata(handles.view_ax,'event_handle');
subject=getappdata(handles.view_ax,'subject');
if ~isempty(event_handle) && ~isempty(event) && event~=hObject.UserData
 % set to normal again
 set(event_handle,'Color','r','LineWidth',0.5)
end
% unset marker
marker=getappdata(handles.view_ax,'marker');
marker_handle=getappdata(handles.view_ax,'marker_handle');
if ~isempty(marker_handle)
 delete(marker_handle)
 setappdata(handles.view_ax,'marker_handle',[]);
end
setappdata(handles.view_ax,'marker',[]);
setappdata(handles.view_ax,'event',hObject.UserData);
event=hObject.UserData;
setappdata(handles.view_ax,'event_handle',hObject);
set(handles.marking_remove,'Enable','On')
update_event_controls(handles)
% populate event name dropdown
evt_IDs={'Spike 1';'Spike 2';'Spike 3'};
if isfield(subject,'evt_IDs')
 evt_IDs=sort(unique(vertcat(evt_IDs,subject.evt_IDs)));
end
% add fixed at the beginning or end
if ~any(strcmpi('',evt_IDs))
 evt_IDs=vertcat({''},evt_IDs);
end
if ~any(strcmpi('Unclear',evt_IDs))
 evt_IDs=vertcat(evt_IDs,{'Unclear'});
end
if ~any(strcmpi('...custom',evt_IDs))
 evt_IDs=vertcat(evt_IDs,{'...custom'});
end
% add IDs
set(handles.event_name,'String',evt_IDs)
% pick the right one
if isfield(subject,'evt_IDs')
 selected=find(strcmpi(subject.evt_IDs{event},evt_IDs));
 if isempty(selected) && ~isempty(subject.evt_IDs{event})
  % seems custom
  set(handles.event_name_custom,'String',subject.evt_IDs{event})
  set(handles.event_name_custom,'Enable','On')
  selected=length(evt_IDs); % always lasr
 end
 set(handles.event_name,'Value',selected)
end
set(handles.event_name,'Enable','On')
% call event redraw
update_event(handles)


% --- Executes on selection change in event_name.
function event_name_Callback(hObject, eventdata, handles)
% hObject    handle to event_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns event_name contents as cell array
%        contents{get(hObject,'Value')} returns selected item from event_name

event=getappdata(handles.view_ax,'event');
event_handle=getappdata(handles.view_ax,'event_handle');
subject=getappdata(handles.view_ax,'subject');
contents = cellstr(get(hObject,'String'));
selected=contents{get(hObject,'Value')};
if ~isfield(subject,'evt_IDs')
 subject.evt_IDs=cell(length(subject.evt_timings_seconds),1);
end
if strcmpi('...custom',selected)
 % enable custom name
 set(handles.event_name_custom,'Enable','On')
else
 set(handles.event_name_custom,'Enable','Off')
 if ~isempty(event) && ~strcmpi(subject.evt_IDs{event},selected)
  % event mode - seems changed
  subject.evt_IDs{event}=selected;
  setappdata(handles.view_ax,'subject',subject);
 end
end
update_view(handles)


% --- Executes during object creation, after setting all properties.
function event_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function event_name_custom_Callback(hObject, eventdata, handles)
% hObject    handle to event_name_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of event_name_custom as text
%        str2double(get(hObject,'String')) returns contents of event_name_custom as a double

event=getappdata(handles.view_ax,'event');
subject=getappdata(handles.view_ax,'subject');
if ~isempty(event)
 % event mode - update ID and view
 subject=getappdata(handles.view_ax,'subject');
 subject.evt_IDs{event}=get(hObject,'String');
 setappdata(handles.view_ax,'subject',subject);
 update_view(handles)
end

% --- Executes during object creation, after setting all properties.
function event_name_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to event_name_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in vertical_align.
function vertical_align_Callback(hObject, eventdata, handles)
% hObject    handle to vertical_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vertical_align contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vertical_align

% just call update - will take care of current state
update_view(handles)


% --- Executes during object creation, after setting all properties.
function vertical_align_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vertical_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

data=getappdata(handles.view_ax,'data');
subject=getappdata(handles.view_ax,'subject');
active_trial=getappdata(handles.view_ax,'active_trial');
trial_substep=getappdata(handles.view_ax,'trial_substep');
xscale=getappdata(handles.view_ax,'xscale');

if eventdata.VerticalScrollCount>0
 % go forward
 for i=1:eventdata.VerticalScrollCount
  trial_substep=trial_substep+1;
  if trial_substep>1/xscale
   active_trial=active_trial+1;
   trial_substep=1;
  end
  if active_trial>length(data.trial)
   active_trial=length(data.trial);
   trial_substep=1/xscale;
  end
 end
end
if eventdata.VerticalScrollCount<0
 % go back
 for i=1:-eventdata.VerticalScrollCount
  trial_substep=trial_substep-1;
  if trial_substep<=0
   active_trial=active_trial-1;
   trial_substep=1/xscale;
  end
  if active_trial<1
   active_trial=1;
   trial_substep=1;
  end
 end
end
setappdata(handles.view_ax,'active_trial',active_trial);
setappdata(handles.view_ax,'trial_substep',trial_substep);
setappdata(handles.view_ax,'event',[]);
update_view(handles)


% --- Executes on mouse press over axes background.
function overview_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to overview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% this is for click on the overview marker
setappdata(handles.view_ax,'active_trial',floor(eventdata.IntersectionPoint(1)));
update_view(handles)


% --- Executes on button press in event_forward.
function event_forward_Callback(hObject, eventdata, handles)
% hObject    handle to event_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
event=getappdata(handles.view_ax,'event');
subject=getappdata(handles.view_ax,'subject');
if isfield(subject,'evt_timings_seconds') && ~isempty(event)
 [evt_s, evt_i]=sort(subject.evt_timings_seconds);
 this_i=find(evt_i==event);
 if this_i<length(evt_s)
  event=evt_i(this_i+1);
 end
 setappdata(handles.view_ax,'event',event);
 setappdata(handles.view_ax,'event_handle',[]);
 update_event(handles)
 update_event_controls(handles)
 set_trial_focus(handles,subject.evt_timings_seconds(event))
end

% --- Executes on button press in event_backward.
function event_backward_Callback(hObject, eventdata, handles)
% hObject    handle to event_backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
event=getappdata(handles.view_ax,'event');
subject=getappdata(handles.view_ax,'subject');
if isfield(subject,'evt_timings_seconds') && ~isempty(event)
 [evt_s, evt_i]=sort(subject.evt_timings_seconds);
 this_i=find(evt_i==event);
 if this_i>1
  event=evt_i(this_i-1);
 end
 setappdata(handles.view_ax,'event',event);
 setappdata(handles.view_ax,'event_handle',[]);
 update_event(handles)
 update_event_controls(handles)
 set_trial_focus(handles,subject.evt_timings_seconds(event))
end

% --- Executes on button press in event_first.
function event_first_Callback(hObject, eventdata, handles)
% hObject    handle to event_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subject=getappdata(handles.view_ax,'subject');
if isfield(subject,'evt_timings_seconds')
 [~,event]=min(subject.evt_timings_seconds);
 setappdata(handles.view_ax,'event',event);
 setappdata(handles.view_ax,'event_handle',[]);
 update_event(handles)
 update_event_controls(handles)
 set_trial_focus(handles,subject.evt_timings_seconds(event))
end

% --- Executes on button press in event_last.
function event_last_Callback(hObject, eventdata, handles)
% hObject    handle to event_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subject=getappdata(handles.view_ax,'subject');
if isfield(subject,'evt_timings_seconds')
 [~,event]=max(subject.evt_timings_seconds);
 setappdata(handles.view_ax,'event',event);
 setappdata(handles.view_ax,'event_handle',[]);
 update_event(handles)
 update_event_controls(handles)
 set_trial_focus(handles,subject.evt_timings_seconds(event))
end


function set_trial_focus(handles,focus)
% will set the active trial / subtrial on focus
data=getappdata(handles.view_ax,'data');
xscale=getappdata(handles.view_ax,'xscale');
% start be determining trial
active_trial=[];
trial_substep=1;
for i=1:length(data.trial)
 if focus>=data.time{i}(1,1) && focus<=data.time{i}(1,end)
  active_trial=i;
 end
end

if ~isempty(active_trial)
 % check substep
 trial_substep=ceil(nearest(data.time{active_trial},focus)/(length(data.time{active_trial})*xscale));
 setappdata(handles.view_ax,'active_trial',active_trial);
 setappdata(handles.view_ax,'trial_substep',trial_substep);
end
% redraw
update_view(handles);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
win_pos=get(hObject,'Position');
% view panel - set to top
old_pos=get(handles.panel_view,'Position');
set(handles.panel_view,'Position',[old_pos(1),win_pos(4)-46,old_pos(3),old_pos(4)]);
% channel panel - set to top
old_pos=get(handles.panel_channels,'Position');
set(handles.panel_channels,'Position',[old_pos(1),win_pos(4)-114,old_pos(3),old_pos(4)]);
% info panel - set to top, same offset
old_pos=get(handles.panel_info,'Position');
set(handles.panel_info,'Position',[old_pos(1),win_pos(4)-46,old_pos(3),old_pos(4)]);
% ratings panel - set to top, same offset
old_pos=get(handles.panel_ratings,'Position');
set(handles.panel_ratings,'Position',[old_pos(1),win_pos(4)-114,old_pos(3),old_pos(4)]);
% events panel - set to top, same offset
old_pos=get(handles.panel_events,'Position');
set(handles.panel_events,'Position',[old_pos(1),win_pos(4)-114,old_pos(3),old_pos(4)]);
% stimuli panel - set to top, same offset
old_pos=get(handles.panel_stimuli,'Position');
set(handles.panel_stimuli,'Position',[old_pos(1),win_pos(4)-114,old_pos(3),old_pos(4)]);
% main view
old_pos=get(handles.view_ax,'Position');
% check if the panels are off
if strcmpi(get(handles.panel_ratings,'Visible'),'Off') && strcmpi(get(handles.panel_events,'Visible'),'Off')
 % move main view up and increase
 set(handles.view_ax,'Position',[old_pos(1),old_pos(2),win_pos(3)-old_pos(1)-12,win_pos(4)-128]);
else
 set(handles.view_ax,'Position',[old_pos(1),old_pos(2),win_pos(3)-old_pos(1)-12,win_pos(4)-197]);
end
% set tick length - otherwise Matlab will make them bigger by default 
ticl=10/(win_pos(3)+win_pos(4));
set(handles.view_ax,'TickLength',[ticl ticl]);
% main view
old_pos=get(handles.overview,'Position');
set(handles.overview,'Position',[old_pos(1),old_pos(2),win_pos(3)-20,old_pos(4)]);
%deal with event window
if isappdata(handles.view_ax,'evFig')
 evFig=getappdata(handles.view_ax,'evFig');
 evLb=getappdata(handles.view_ax,'evLb');
 win_pos=get(hObject,'Position');
 if ishandle(evFig)
  set(evFig,'Position',[win_pos(1)+win_pos(3)+10,win_pos(2),300,win_pos(4)]);  
  win_pos=get(evFig,'Position');
  if ishandle(evLb)
  set(evLb,'Position',[10 10 win_pos(3)-20  win_pos(4)-20 ]);
  end
 end
end
%deal with layout window
if isappdata(handles.view_ax,'layFig')
 layFig=getappdata(handles.view_ax,'layFig');
 win_pos=get(hObject,'Position');
 if ishandle(layFig)
  this_pos=get(layFig,'Position');
  set(layFig,'Position',[win_pos(1)-10-this_pos(3),win_pos(2)+win_pos(4)-this_pos(4),this_pos(3),this_pos(4)]);  
 end
end



% --- Executes on button press in skip_next_vigilance.
function skip_next_vigilance_Callback(hObject, eventdata, handles)
% hObject    handle to skip_next_vigilance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skip_next_vigilance


function evLb_selection_Callback(handles,src,event)
% callback for Event listbox
if isappdata(handles.view_ax,'evJumplist')
 jumplist=getappdata(handles.view_ax,'evJumplist');
 if ~isempty(jumplist) 
  sel=get(src,'Value');
  % check if in data
  if ~isnan(jumplist(sel,3))
   % check if evt
   if ~isnan(jumplist(sel,2))
    setappdata(handles.view_ax,'event',jumplist(sel,2));
    setappdata(handles.view_ax,'event_handle',[]);
    update_event(handles)
    update_event_controls(handles)
   end
   set_trial_focus(handles,jumplist(sel,1));
  end
 end
end


% --- Executes on selection change in stimuli.
function stimuli_Callback(hObject, eventdata, handles)
% hObject    handle to stimuli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stimuli contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stimuli

data=getappdata(handles.view_ax,'data');
subject=getappdata(handles.view_ax,'subject');
active_trial=getappdata(handles.view_ax,'active_trial');
switch get(hObject,'Value')
 case 1
  data.trial_markings{active_trial,4}='eyes_closed';
 case 2
  data.trial_markings{active_trial,4}='eyes_open';
 case 3
  data.trial_markings{active_trial,4}='HV';
 case 4
  data.trial_markings{active_trial,4}='postHV';
 case 5
  data.trial_markings{active_trial,4}='PS';
 case 6
  data.trial_markings{active_trial,4}='';
end
setappdata(handles.view_ax,'data',data);
% update bg_cols
bg_cols=make_bg_cols(data,subject); 
setappdata(handles.view_ax,'bg_cols',bg_cols);
update_view(handles)


% --- Executes during object creation, after setting all properties.
function stimuli_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimuli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in skip_next_stimuli.
function skip_next_stimuli_Callback(hObject, eventdata, handles)
% hObject    handle to skip_next_stimuli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skip_next_stimuli


% --- Executes on button press in hide_bad_channels.
function hide_bad_channels_Callback(hObject, eventdata, handles)
% hObject    handle to hide_bad_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hide_bad_channels

 select_channels(handles)
 update_view(handles)


% --- Executes on button press in toggle_channel.
function toggle_channel_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
active_channel=getappdata(handles.view_ax,'active_channel'); 
if ~isempty(active_channel)
 data=getappdata(handles.view_ax,'data');
 if any(strcmp(active_channel,data.bad_channels))
  %remove from bad
  idx=find(strcmp(active_channel,data.bad_channels));
  data.bad_channels=[data.bad_channels(1:idx-1) data.bad_channels(idx+1:end)];
  fprintf('Flagging channel=%s as GOOD\n',active_channel)
 else
  data.bad_channels{end+1}=active_channel;
  fprintf('Flagging channel=%s as BAD\n',active_channel)
 end
 setappdata(handles.view_ax,'data',data);
 select_channels(handles)
 update_view(handles)
 % update layout
 if isappdata(handles.view_ax,'layFig')
  layFig=getappdata(handles.view_ax,'layFig');
  if ishandle(layFig)
   update_layout_channels(handles)
  end
 end
end


% --- Executes on button press in dash_bad_channels.
function dash_bad_channels_Callback(hObject, eventdata, handles)
% hObject    handle to dash_bad_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dash_bad_channels

 update_view(handles)

 
 
 
% draw layout 
function draw_layout(handles)
subject=getappdata(handles.view_ax,'subject');
active_channel=getappdata(handles.view_ax,'active_channel');
% udpate event window, if any
if isappdata(handles.view_ax,'layFig')
 layFig=getappdata(handles.view_ax,'layFig');
 % set colors
 data=getappdata(handles.view_ax,'data');
 elec_col=repmat({[0 0.5 0]},length(subject.layout.label),1);
 elec_bold=repmat({'normal'},length(subject.layout.label),1);
 if isfield(data,'bad_channels')
  for i=1:length(subject.layout.label)
   if any(strcmp(subject.layout.label{i},data.bad_channels))
    elec_col{i}='r';
   end
  end
 end
 if ~isempty(active_channel)
  elec_bold{find(strcmp(subject.layout.label,active_channel),1)}='bold';
 end
 
 % plot
 figure(layFig)
 clf(layFig)
 nmri_plot_layout(subject.layout,'box','no','point','no','channel_hit_callback',@(hObject,eventdata)eeg_score('channel_layout_click_Callback',hObject,eventdata,handles))
 zoom(layFig,1.1)
 update_layout_channels(handles)
 figure(handles.figure1)
end 


% reformat the channels in layout
function update_layout_channels(handles)
subject=getappdata(handles.view_ax,'subject');
active_channel=getappdata(handles.view_ax,'active_channel');
data=getappdata(handles.view_ax,'data');
if isappdata(handles.view_ax,'layFig')
 layFig=getappdata(handles.view_ax,'layFig');
 
 items=findobj(layFig,'Type','Text');
  % make colors and bold
 data=getappdata(handles.view_ax,'data');
 elec_col=repmat({[0 0.5 0]},length(subject.layout.label),1);
 elec_bold=repmat({'normal'},length(subject.layout.label),1);
 if isfield(data,'bad_channels')
  for i=1:length(subject.layout.label)
   if any(strcmp(subject.layout.label{i},data.bad_channels))
    elec_col{i}='r';
   elseif ~any(strcmpi(subject.layout.label{i},data.label))
    % check if in data
    elec_col{i}=[0.8 0.8 0.8]; % very grey
   end
  end
 end
 if ~isempty(active_channel)
  elec_bold{find(strcmp(subject.layout.label,active_channel),1)}='bold';
 end
 
 % now loop all text items
 for i=1:numel(items)
  idx=find(strcmp(items(i).String,subject.layout.label));
  if ~isempty(idx)
   set(items(i),'Color',elec_col{idx})
   set(items(i),'FontWeight',elec_bold{idx})
  else
   set(items(i),'Color','b')
   set(items(i),'FontWeight','normal')
  end
  
 end
 
end


% --- Executes on button press in layout
function channel_layout_click_Callback(hObject, eventdata, handles)
active_channel=getappdata(handles.view_ax,'active_channel');
data=getappdata(handles.view_ax,'data');

fprintf('Clicked on channel=%s\n',hObject.String)
if strcmp(hObject.String,active_channel)
 active_channel=[];
 set(handles.toggle_channel,'Enable','Off')
else
 if any(strcmp(hObject.String,data.label))
  active_channel=hObject.String;
  set(handles.toggle_channel,'Enable','On')
 end
end
setappdata(handles.view_ax,'active_channel',active_channel); 
update_view(handles)
% update layout
if isappdata(handles.view_ax,'layFig')
 layFig=getappdata(handles.view_ax,'layFig');
 if ishandle(layFig)
  update_layout_channels(handles)
 end
end


 


% --- Executes on button press in layout_show.
function layout_show_Callback(hObject, eventdata, handles)
% hObject    handle to layout_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subject=getappdata(handles.view_ax,'subject'); 

openLay=true;
if isappdata(handles.view_ax,'layFig')
 layFig=getappdata(handles.view_ax,'layFig');
 if ishandle(layFig)
  %close
  close(layFig)
  openLay=false;
 end
end
if openLay
 % instead raise the layout
 % make the layout
 if ~isstruct(subject.layout)
  lcfg=[];
  lcfg.layout=subject.layout;
  subject.layout = ft_prepare_layout(lcfg);
 end
 
 % remove scale and comment
 remove_ch=strcmpi(subject.layout.label,'scale')|strcmpi(subject.layout.label,'comnt');
 if ~isempty(remove_ch)
  subject.layout.label= subject.layout.label(~remove_ch,:);
  subject.layout.pos= subject.layout.pos(~remove_ch,:);
  subject.layout.width= subject.layout.width(~remove_ch,:);
  subject.layout.height= subject.layout.height(~remove_ch,:);  
 end
 
 win_pos=get(handles.figure1,'Position');
 if numel(subject.layout.label)>32
  sizeL=600;
 else
  sizeL=300;
 end
 
 layFig=figure('Toolbar','none','DockControls','Off','MenuBar','none','Name','Layout','NumberTitle','off',...
  'Position',[win_pos(1)-12-sizeL,win_pos(2)+win_pos(4)-sizeL,sizeL,sizeL]);
 setappdata(handles.view_ax,'layFig',layFig);
 setappdata(handles.view_ax,'subject',subject);
 
 draw_layout(handles)
 
 %give back the focus to hObject
 set(groot,'CurrentFigure',handles.figure1);
 % make the main figure current again
 figure(handles.figure1)
end

% --- Executes on button press in event_list_show.
function event_list_show_Callback(hObject, eventdata, handles)
% hObject    handle to event_list_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subject=getappdata(handles.view_ax,'subject'); 

openList=true;
if isappdata(handles.view_ax,'evFig')
 evFig=getappdata(handles.view_ax,'evFig');
 if ishandle(evFig)
  %close
  close(evFig)
  openList=false;
 end
end
if openList
 % instead raise the Event Lsit
 win_pos=get(handles.figure1,'Position');
 evFig=figure('Toolbar','none','DockControls','Off','MenuBar','none','Name','Events/Infos','NumberTitle','off',...
  'Position',[win_pos(1)+win_pos(3)+16,win_pos(2),300,win_pos(4)]);
 %give back the focus to hObject
 set(groot,'CurrentFigure',handles.figure1);
 % add a list
 win_pos=get(evFig,'Position');
 evLb=uicontrol(evFig,'Style','listbox','Position',[10 10 win_pos(3)-20  win_pos(4)-20 ],'Callback',@(src,event)evLb_selection_Callback(handles,src,event));
 setappdata(handles.view_ax,'evFig',evFig);
 setappdata(handles.view_ax,'evLb',evLb);
 % make the main figure current again
 figure(handles.figure1)
 update_evt_list(handles)
end
