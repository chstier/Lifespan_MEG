function varargout = ica_selection(varargin)
% ICA_SELECTION MATLAB code for ica_selection.fig
%      ICA_SELECTION, by itself, creates a new ICA_SELECTION or raises the existing
%      singleton*.
%
%      H = ICA_SELECTION returns the handle to a new ICA_SELECTION or the handle to
%      the existing singleton*.
%
%      ICA_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ICA_SELECTION.M with the given input arguments.
%
%      ICA_SELECTION('Property','Value',...) creates a new ICA_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ica_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ica_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ica_selection

% Last Modified by GUIDE v2.5 14-Mar-2017 13:40:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ica_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @ica_selection_OutputFcn, ...
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


% --- Executes just before ica_selection is made visible.
function ica_selection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ica_selection (see VARARGIN)

% Choose default command line output for ica_selection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% check if we have subject struct
if length(varargin)>1 && ~isempty(varargin{1}) && ~isempty(varargin{2}) 
 comp=varargin{2};
 cfg=varargin{1};
else
 error('Need component struct from Fieldtrip')
end
if ~isfield(comp,'topo') || ~isfield(comp,'time') || ~isfield(comp,'trial')
 error('Component struct seems incomplete / not valid')
end

if ~isfield(cfg,'layout') 
 error('Layout is a required field for cfg')
end


% make layout
if ~isstruct(cfg.layout)
 cfg.layout = ft_prepare_layout(cfg);
end
if ~isfield(cfg,'selected') || length(cfg.selected)~=size(comp.topo,2)
 % if not given, set to all ones
 cfg.selected=ones(size(comp.topo,2),1);
end

setappdata(handles.comp_1,'comp',comp)
setappdata(handles.comp_1,'cfg',cfg)
setappdata(handles.comp_1,'start_comp',1)
setappdata(handles.comp_1,'active_trial',floor(length(comp.trial)/2))
setappdata(handles.comp_1,'comp_toggle',cfg.selected)

% setup sliders
set(handles.trial_slider,'Max',length(comp.trial))
set(handles.trial_slider,'Min',1)
set(handles.trial_slider,'SliderStep',[1/(length(comp.trial)-1) 1/(length(comp.trial)-1)])
set(handles.trial_slider,'Value',1)

set(handles.comp_slider,'Max',size(comp.topo,2)-4)
set(handles.comp_slider,'Min',1)
set(handles.comp_slider,'SliderStep',[(5/(size(comp.topo,2)-5)) (5/(size(comp.topo,2)-5))])
set(handles.comp_slider,'Value',size(comp.topo,2)-4)

update_view(handles)
% UIWAIT makes ica_selection wait for user response (see UIRESUME)
uiwait(handles.ica_selection_fig);


function update_view(handles)
comp=getappdata(handles.comp_1,'comp');
cfg=getappdata(handles.comp_1,'cfg');
start_comp=getappdata(handles.comp_1,'start_comp');
active_trial=getappdata(handles.comp_1,'active_trial');
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
[seldat, sellay] = match_str(comp.topolabel, cfg.layout.label); 
if isempty(seldat)
 error('labels in data and labels in layout do not match');
end
chanX = cfg.layout.pos(sellay,1);
chanY = cfg.layout.pos(sellay,2);
chanLabels = cfg.layout.label(sellay);
opt = {'outline',cfg.layout.outline,...
'interpmethod', 'v4',...
'interplim','mask',...
'mask',cfg.layout.mask};

% cycle through components
for n=1:5
 comp_n=n+start_comp-1;
 set(handles.(['panel_' num2str(n)]),'Title',['Component #' num2str(comp_n)])
 dat = comp.topo(seldat,comp_n);
 axes(handles.(['comp_' num2str(n)]));
 cla
 ft_plot_topo(chanX,chanY,dat,opt{:});
 % newer Fieldtrips have a function call change
 if exist('ft_plot_layout','file')==2
  ft_plot_layout(cfg.layout,'box','no','label','no','point','no')
 else
  ft_plot_lay(cfg.layout,'box','no','label','no','point','no')
 end
 set(handles.(['comp_' num2str(n)]),'HitTest','Off')
 set(allchild(handles.(['comp_' num2str(n)])),'HitTest','Off')
 % timecourse
 axes(handles.(['timecourse_' num2str(n)]));
 cla
 plot(comp.time{active_trial}(1:end),comp.trial{active_trial}(comp_n,:));
 yscale=range(comp.trial{active_trial}(comp_n,:));
 ylim([-yscale yscale])
 set(handles.(['timecourse_' num2str(n)]),'HitTest','Off','YTickLabel','')
 set(allchild(handles.(['timecourse_' num2str(n)])),'HitTest','Off')
 % comp toggel
 if ~comp_toggle(comp_n,1)
  % reject
  set(handles.(['panel_' num2str(n)]),'BackgroundColor',[1 0 0]);
  set(handles.(['keep_' num2str(n)]),'Value',2);
 else
  % keep
  set(handles.(['panel_' num2str(n)]),'BackgroundColor',[0.94 0.94 0.94]);
  set(handles.(['keep_' num2str(n)]),'Value',1);
 end
end
set(handles.trial_slider,'Value',active_trial)
set(handles.comp_slider,'Value',size(comp.topo,2)-start_comp-3)

% --- Outputs from this function are returned to the command line.
function varargout = ica_selection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Return selected components
varargout{1} = getappdata(handles.comp_1,'comp_toggle');
delete(hObject);

% --- Executes when user attempts to close ica_selection_fig.
function ica_selection_fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to ica_selection_fig (see GCBO)
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


% --- Executes on button press in keep_1.
function keep_1_Callback(hObject, eventdata, handles)
% hObject    handle to keep_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(hObject,'Value');
if current==2
 set(handles.panel_1,'BackgroundColor',[1 0 0]);
else
 set(handles.panel_1,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp,1)=current==1;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)


% --- Executes on button press in keep_5.
function keep_5_Callback(hObject, eventdata, handles)
% hObject    handle to keep_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(hObject,'Value');
if current==2
 set(handles.panel_5,'BackgroundColor',[1 0 0]);
else
 set(handles.panel_5,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+4,1)=current==1;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)

% --- Executes on button press in keep_4.
function keep_4_Callback(hObject, eventdata, handles)
% hObject    handle to keep_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(hObject,'Value');
if current==2
 set(handles.panel_4,'BackgroundColor',[1 0 0]);
else
 set(handles.panel_4,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+3,1)=current==1;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)

% --- Executes on button press in keep_3.
function keep_3_Callback(hObject, eventdata, handles)
% hObject    handle to keep_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(hObject,'Value');
if current==2
 set(handles.panel_3,'BackgroundColor',[1 0 0]);
else
 set(handles.panel_3,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+2,1)=current==1;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)

% --- Executes on button press in keep_2.
function keep_2_Callback(hObject, eventdata, handles)
% hObject    handle to keep_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(hObject,'Value');
if current==2
 set(handles.panel_2,'BackgroundColor',[1 0 0]);
else
 set(handles.panel_2,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+1,1)=current==1;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)

% --- Executes on key press with focus on ica_selection_fig or any of its controls.
function ica_selection_fig_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to ica_selection_fig (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
comp=getappdata(handles.comp_1,'comp');
cfg=getappdata(handles.comp_1,'cfg');
start_comp=getappdata(handles.comp_1,'start_comp');
active_trial=getappdata(handles.comp_1,'active_trial');

% master key switch
switch(eventdata.Key)
 case 'rightarrow'
  % go trial forward
  active_trial=active_trial+1;
  if active_trial>length(comp.trial)
   active_trial=length(comp.trial);
  end
  setappdata(handles.comp_1,'active_trial',active_trial);
  update_view(handles)
 case 'leftarrow'
  % go trial back
  active_trial=active_trial-1;
  if active_trial<1
   active_trial=1;
  end
  setappdata(handles.comp_1,'active_trial',active_trial);
  update_view(handles)
 case 'uparrow'
  % go components down
  start_comp=start_comp-5;
  if start_comp<1
   start_comp=1;
  end
  setappdata(handles.comp_1,'start_comp',start_comp);
  update_view(handles)
 case 'downarrow'
  % go components up
  start_comp=start_comp+5;
  if start_comp>(size(comp.topo,2)-4)
   start_comp=(size(comp.topo,2)-4);
  end
  setappdata(handles.comp_1,'start_comp',start_comp);
  update_view(handles)
end


% --------------------------------------------------------------------
function panel_1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to panel_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(handles.keep_1,'Value');
set(handles.keep_1,'Value',3-current);
if current==1
 set(hObject,'BackgroundColor',[1 0 0]);
else
 set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp,1)=current==2;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)


% --------------------------------------------------------------------
function panel_2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to panel_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(handles.keep_2,'Value');
set(handles.keep_2,'Value',3-current);
if current==1
 set(hObject,'BackgroundColor',[1 0 0]);
else
 set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+1,1)=current==2;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)


% --------------------------------------------------------------------
function panel_3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to panel_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(handles.keep_3,'Value');
set(handles.keep_3,'Value',3-current);
if current==1
 set(hObject,'BackgroundColor',[1 0 0]);
else
 set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+2,1)=current==2;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)


% --------------------------------------------------------------------
function panel_4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to panel_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(handles.keep_4,'Value');
set(handles.keep_4,'Value',3-current);
if current==1
 set(hObject,'BackgroundColor',[1 0 0]);
else
 set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+3,1)=current==2;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)


% --------------------------------------------------------------------
function panel_5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to panel_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current=get(handles.keep_5,'Value');
set(handles.keep_5,'Value',3-current);
if current==1
 set(hObject,'BackgroundColor',[1 0 0]);
else
 set(hObject,'BackgroundColor',[0.94 0.94 0.94]);
end
comp_toggle=getappdata(handles.comp_1,'comp_toggle');
start_comp=getappdata(handles.comp_1,'start_comp');
comp_toggle(start_comp+4,1)=current==2;
setappdata(handles.comp_1,'comp_toggle',comp_toggle)


% --- Executes on slider movement.
function comp_slider_Callback(hObject, eventdata, handles)
% hObject    handle to comp_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value=round(get(hObject,'Value'));
max=round(get(hObject,'Max'));
setappdata(handles.comp_1,'start_comp',max-value+1);
update_view(handles)

% --- Executes during object creation, after setting all properties.
function comp_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comp_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function trial_slider_Callback(hObject, eventdata, handles)
% hObject    handle to trial_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
active_trial=round(get(hObject,'Value'));
setappdata(handles.comp_1,'active_trial',active_trial);
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


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ica_selection_fig_CloseRequestFcn(handles.ica_selection_fig, [], handles)
