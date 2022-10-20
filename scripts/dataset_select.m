function varargout = dataset_select(varargin)
% DATASET_SELECT MATLAB code for dataset_select.fig
%      DATASET_SELECT, by itself, creates a new DATASET_SELECT or raises the existing
%      singleton*.
%
%      H = DATASET_SELECT returns the handle to a new DATASET_SELECT or the handle to
%      the existing singleton*.
%
%      DATASET_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATASET_SELECT.M with the given input arguments.
%
%      DATASET_SELECT('Property','Value',...) creates a new DATASET_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dataset_select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dataset_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dataset_select

% Last Modified by GUIDE v2.5 16-Mar-2017 14:16:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dataset_select_OpeningFcn, ...
                   'gui_OutputFcn',  @dataset_select_OutputFcn, ...
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


% --- Executes just before dataset_select is made visible.
function dataset_select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dataset_select (see VARARGIN)

% Choose default command line output for dataset_select
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% check if we have subject struct
if length(varargin)>0 && ~isempty(varargin{1})
 subject=varargin{1};
else
 error('No subject info given, need subject struct')
end
if ~isfield(subject,'id') || ~isfield(subject,'exam_id') || ~isfield(subject,'analysis_dir')
 error('Subject struct seems incomplete / not valid')
end
txt={['ID=' subject.id ', Exam_ID=' subject.exam_id],['Analyis_Dir=' subject.analysis_dir]};
set(handles.subject_info,'String',txt)

% now check all datasets (potentially) available
pot_datasets={'dws_filt_dataset','clean_dataset','cleanICA_dataset'};
pot_datasets_labels={'Downsampled/Filtered','Cleaned','Cleaned+ICA'};
datasets={};
datasets_labels={};
for d=1:length(pot_datasets)
 if (isfield(subject,pot_datasets{d}) && exist(subject.(pot_datasets{d}),'file')) 
  datasets{end+1}=subject.(pot_datasets{d});
  datasets_labels(end+1)=pot_datasets_labels(d);  
 end
end


% now make data handing structure
all_data=cell(length(datasets),1);
for d=1:length(datasets)
 % load subject info for all now
 all_data{d}=load(datasets{d},'subject');
 % check for re-root
 if (~strcmp(datasets{d}(1:length(all_data{d}.subject.analysis_dir)),all_data{d}.subject.analysis_dir))
  % there has been a move, now try to determine the true analysis dir
  analysis_dir=dirname(dirname(dirname(datasets{d})));
  all_data{d}.subject=nmri_reroot_subject(all_data{d}.subject,analysis_dir); 
 end
 all_data{d}.label=datasets_labels{d};
 all_data{d}.dataset=datasets{d};
end
setappdata(handles.datasets_overview,'all_data',all_data)
setappdata(handles.datasets_overview,'selected_dataset',[])

% fill list initially
update_list(handles)

% UIWAIT makes dataset_select wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function update_list(handles)
all_data=getappdata(handles.datasets_overview,'all_data');
selected_dataset=getappdata(handles.datasets_overview,'selected_dataset');
% populate table
dataset_list={};
for d=1:length(all_data)
 goodTrials=[];
 goodChannels=[];
 subject=all_data{d}.subject;
 if isfield(subject,'evt_timings_seconds')
  evts=length(subject.evt_timings_seconds);
 else
  evts='none';
 end
 if isfield(subject,'data_info')
  % if we have this (saved by nmri_write_dataset) we take this, much
  % faster and smaller
  if isfield(subject.data_info,'trial_markings')
   marks=sprintf('Wake: %d, Sleep: %d, Bad: %d, ToDo %d',sum(strcmp(subject.data_info.trial_markings(:,1),'wake')),(sum(strcmp(subject.data_info.trial_markings(:,1),'sleep1'))...
+sum(strcmp(subject.data_info.trial_markings(:,1),'sleep2'))+sum(strcmp(subject.data_info.trial_markings(:,1),'sleep3'))+sum(strcmp(subject.data_info.trial_markings(:,1),'rem'))),...
sum(cellfun(@(x) islogical(x) && x==false,subject.data_info.trial_markings(:,2))),sum(cellfun(@(x) isempty(x),subject.data_info.trial_markings(:,1))));
  else
   marks='none';
  end
  if isfield(subject.data_info,'nchannels')
   nchan=subject.data_info.nchannels;
  else
   nchan='N/A';
  end
  if isfield(subject.data_info,'ntrials')
   ntrials=subject.data_info.ntrials;
  else
   ntrials='N/A';
  end
  if isfield(subject.data_info,'fsample') && isfield(subject.data_info,'trial_markings_sampleinfo')
   % fake data.time
   subject.data_info.time=subject.data_info.trial_markings_sampleinfo(:,2);
   [ goodTrials, ~ ] = nmri_trial_selector(subject,subject.data_info,[],0);
  else
   goodTrials=[];
  end
  
  if isfield(subject.data_info,'bad_channels') && isfield(subject.data_info,'nchannels')
   goodChannels=subject.data_info.nchannels-length(subject.data_info.bad_channels);
  end
  
 else
  % then we need the real data -- slow, only do on demand
  marks='[click to check data]';
  nchan='[click]';
  ntrials='[click]';
 end
 % modified
 if isfield(all_data{d},'modified')
  modified=datestr(all_data{d}.modified,'HH:MM');
 else
  modified='';
 end
 
  % goodtrials
 if ~isempty(goodTrials)
  goodT=num2str(length(goodTrials));
 else
  goodT='N/A';
 end
 
 %good channels
 if ~isempty(goodChannels)
  goodCh=num2str(goodChannels);
 else
  goodCh='N/A';
 end
 
 
 % now make info fields
 dataset_list=[dataset_list; {all_data{d}.label,nchan,ntrials,evts,marks,modified,goodT,goodCh}];
 
 % now check controls
 if ~isempty(selected_dataset) && (selected_dataset==d)
  if isfield(all_data{d},'modified')
   set(handles.save,'Enable','On')
   set(handles.save_as,'Enable','On')
  else
   set(handles.save,'Enable','Off')
   set(handles.save_as,'Enable','Off')   
  end
  set(handles.view,'Enable','On')
  set(handles.delete,'Enable','On')
 end
end
if isempty(dataset_list)
 dataset_list={'No datasets found','','','','',''};
end
set(handles.datasets_overview,'Data',dataset_list)

% --- Outputs from this function are returned to the command line.
function varargout = dataset_select_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in view.
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open the current selection
all_data=getappdata(handles.datasets_overview,'all_data');
selected_dataset=getappdata(handles.datasets_overview,'selected_dataset');
if ~isempty(selected_dataset)
 if ~isfield(all_data{selected_dataset},'data')
  % load data
  load_dataset(handles,selected_dataset)
  all_data=getappdata(handles.datasets_overview,'all_data');
 end
 [data, subject]=eeg_score([],all_data{selected_dataset}.data,all_data{selected_dataset}.subject);
 if ~isequal(data,all_data{selected_dataset}.data) || ~isequal(subject,all_data{selected_dataset}.subject)
  % update if different
  all_data{selected_dataset}.data=data;
  all_data{selected_dataset}.subject=subject;
  all_data{selected_dataset}.modified=datetime;
  setappdata(handles.datasets_overview,'all_data',all_data);
  clear data subject
  update_list(handles);
 end
end


% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
all_data=getappdata(handles.datasets_overview,'all_data');
selected_dataset=getappdata(handles.datasets_overview,'selected_dataset');
button=questdlg({'Do you really want to delete this file?',['Dataset = ' all_data{selected_dataset}.label]},'Confirm Delete','No');
if strcmpi(button,'Yes')
 h=msgbox({['Deleting Dataset = ' all_data{selected_dataset}.label],'Please wait...'});
 delete(all_data{selected_dataset}.dataset);
 all_data(selected_dataset)=[];
 selected_dataset=[];
 setappdata(handles.datasets_overview,'all_data',all_data);
 setappdata(handles.datasets_overview,'selected_dataset',selected_dataset);
 update_list(handles)
 try
  delete(h)
 catch
 end
end


% --- Executes on button press in save_as.
function save_as_Callback(hObject, eventdata, handles)
% hObject    handle to save_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
all_data=getappdata(handles.datasets_overview,'all_data');
selected_dataset=getappdata(handles.datasets_overview,'selected_dataset');
[target_file, target_path]=uiputfile('*.mat');
target=fullfile(target_path,target_file);
if ischar(target_file)
 h=msgbox({['Saving Dataset = ' all_data{selected_dataset}.label],'Please wait...'});
 nmri_write_dataset(target,all_data{selected_dataset}.data,all_data{selected_dataset}.subject)
 all_data{selected_dataset}.saved=datetime;
 if isfield(all_data{selected_dataset},'modified')
  all_data{selected_dataset}=rmfield(all_data{selected_dataset},'modified');
 end
 setappdata(handles.datasets_overview,'all_data',all_data);
 update_list(handles)
 try
  delete(h)
 catch
 end
end


% --- Executes when selected cell(s) is changed in datasets_overview.
function datasets_overview_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to datasets_overview (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
all_data=getappdata(handles.datasets_overview,'all_data');
selected_datasets=unique(eventdata.Indices(:,1));
if length(selected_datasets)==1
 % single dataset
 set(handles.view,'Enable','On')
 set(handles.delete,'Enable','On')
 setappdata(handles.datasets_overview,'selected_dataset',selected_datasets(1))
 % check if we want to load full data
 if eventdata.Indices(1,2)~=1 && ~isfield(all_data{selected_datasets(1)}.subject,'data_info')
  % load data
  load_dataset(handles,selected_datasets(1))
 end 
elseif length(selected_datasets)>1
 set(handles.delete,'Enable','On')
 set(handles.save,'Enable','Off') 
 set(handles.save_as,'Enable','Off') 
 set(handles.view,'Enable','Off') 
else
 %nothing
 set(handles.view,'Enable','Off')
 set(handles.delete,'Enable','Off')
 set(handles.save,'Enable','Off') 
 set(handles.save_as,'Enable','Off') 
end
update_list(handles);

function load_dataset(handles,selected_dataset)
all_data=getappdata(handles.datasets_overview,'all_data');
h=msgbox({['Loading Dataset = ' all_data{selected_dataset}.label],'Please wait...'});
data=load(all_data{selected_dataset}.dataset,'data');
all_data{selected_dataset}.data=data.data;
clear data
all_data{selected_dataset}.subject.data_info.ntrials=size(all_data{selected_dataset}.data.trial,2);
all_data{selected_dataset}.subject.data_info.nchannels=size(all_data{selected_dataset}.data.trial{1},1);
all_data{selected_dataset}.loaded=datetime;
all_data{selected_dataset}.subject.data_info.nchannels=size(all_data{selected_dataset}.data.trial{1},1);
all_data{selected_dataset}.subject.data_info.fsample=all_data{selected_dataset}.data.fsample;
if isfield(all_data{selected_dataset}.data,'bad_channels')
 all_data{selected_dataset}.subject.data_info.bad_channels=all_data{selected_dataset}.data.bad_channels;
end
setappdata(handles.datasets_overview,'all_data',all_data);
try
 delete(h)
catch
end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
all_data=getappdata(handles.datasets_overview,'all_data');
selected_dataset=getappdata(handles.datasets_overview,'selected_dataset');
button=questdlg({'Do you really want to save/overwrite?',['Dataset = ' all_data{selected_dataset}.label]},'Confirm Save','No');
if strcmpi(button,'Yes')
 h=msgbox({['Saving Dataset = ' all_data{selected_dataset}.label],'Please wait...'});
 nmri_write_dataset(all_data{selected_dataset}.dataset,all_data{selected_dataset}.data,all_data{selected_dataset}.subject)
 all_data{selected_dataset}.saved=datetime;
 if isfield(all_data{selected_dataset},'modified')
  all_data{selected_dataset}=rmfield(all_data{selected_dataset},'modified');
 end
 setappdata(handles.datasets_overview,'all_data',all_data);
 update_list(handles)
 try
  delete(h)
 catch
 end
end

% --- Executes during object deletion, before destroying properties.
function view_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function delete_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function save_as_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to save_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function save_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
all_data=getappdata(handles.datasets_overview,'all_data');
% check for unsaved changes
unsaved=0;
for d=1:length(all_data)
 if isfield(all_data{d},'modified')
  unsaved=unsaved+1;
 end
end
if unsaved>0
 button=questdlg({['There are ' num2str(unsaved) ' modified and unsaved datasets.'],'These changes will be lost upon closeing.','Do you really want to close without saving?'},'Confirm Close','No');
 if strcmpi(button,'Yes')
  % alright - so close
  delete(hObject);
 end
else
 % nothing unsaved - close
 delete(hObject);
end 


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
