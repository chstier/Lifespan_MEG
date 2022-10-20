function varargout = mri_view(varargin)
% MRI_VIEW MATLAB code for mri_view.fig
%      MRI_VIEW, by itself, creates a new MRI_VIEW or raises the existing
%      singleton*.
%
%      H = MRI_VIEW returns the handle to a new MRI_VIEW or the handle to
%      the existing singleton*.
%
%      MRI_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRI_VIEW.M with the given input arguments.
%
%      MRI_VIEW('Property','Value',...) creates a new MRI_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mri_view_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mri_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mri_view

% Last Modified by GUIDE v2.5 31-May-2019 16:45:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mri_view_OpeningFcn, ...
                   'gui_OutputFcn',  @mri_view_OutputFcn, ...
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


% --- Executes just before mri_view is made visible.
function mri_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mri_view (see VARARGIN)

% Choose default command line output for mri_view
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% check if we have cfg
if length(varargin)>0 && ~isempty(varargin{1})
 cfg=varargin{1};
else
 cfg=struct('empty',[]);
end
% check if we have mri
if length(varargin)>1
 mri=varargin{2};
else
 mri=[];
end
% check if we have subject
if length(varargin)>2
 subject=varargin{3};
else
 subject=[];
end

% make a focus
if isfield(cfg,'focus')
 focus=cfg.focus;
elseif isfield(mri,'anatomy')
 focus=round(size(mri.anatomy)/2);
end

% read voxel size
P=imatrix(mri.transform);
cfg.vx=P(7:9);

% set the sliders
set(handles.sliderX,'Max',mri.dim(3))
set(handles.sliderX,'Min',1)
set(handles.sliderX,'SliderStep',[1/mri.dim(3) 10/mri.dim(3)])
set(handles.sliderX,'Value',focus(3))

set(handles.sliderY,'Max',mri.dim(2))
set(handles.sliderY,'Min',1)
set(handles.sliderY,'SliderStep',[1/mri.dim(2) 10/mri.dim(2)])
set(handles.sliderY,'Value',focus(2))

set(handles.sliderZ,'Max',mri.dim(1))
set(handles.sliderZ,'Min',1)
set(handles.sliderZ,'SliderStep',[1/mri.dim(1) 10/mri.dim(1)])
set(handles.sliderZ,'Value',focus(1))


% set the window
if ~isfield(cfg,'max') 
 if isfield(mri,'anatomy')
  cfg.max=max(mri.anatomy(:));
 else
  cfg.max=9999;
 end
end
if ~isfield(cfg,'min')
 if isfield(mri,'anatomy')
  cfg.min=min(mri.anatomy(:));
 else
  cfg.min=0;
 end
end

if ~isfield(cfg,'window_min') 
 if isfield(mri,'anatomy')
  cfg.window_min=mean(mri.anatomy(:))-std(mri.anatomy(:));
 else
  cfg.window_min=cfg.min;
 end
end

if ~isfield(cfg,'window_max') 
 if isfield(mri,'anatomy')
  cfg.window_max=mean(mri.anatomy(:))+(2*std(mri.anatomy(:)));
 else
  cfg.window_max=cfg.max;
 end
end

if cfg.window_max>cfg.max
 cfg.window_max=cfg.max;
end
if cfg.window_min<cfg.min
 cfg.window_min=cfg.min;
end


% set the vol tresh

if ~isfield(cfg,'volThresh') 
 cfg.volThresh=mean(mri.anatomy(:));
end

if cfg.volThresh>cfg.max
 cfg.volThresh=cfg.max;
end
if cfg.volThresh<cfg.min
 cfg.volThresh=cfg.min;
end

% set the min-max sliders
set(handles.WindowMin,'Min',cfg.min)
set(handles.WindowMax,'Min',cfg.min)
set(handles.WindowMin,'Max',cfg.max)
set(handles.WindowMax,'Max',cfg.max)

set(handles.WindowMin,'SliderStep',[1/(cfg.max-cfg.min) 10/(cfg.max-cfg.min)])
set(handles.WindowMin,'Value',cfg.window_min)

set(handles.WindowMax,'SliderStep',[1/(cfg.max-cfg.min) 10/(cfg.max-cfg.min)])
set(handles.WindowMax,'Value',cfg.window_max)

set(handles.WindowMaxVal,'String',cfg.window_max)
set(handles.WindowMinVal,'String',cfg.window_min)

% VolThr
set(handles.VolThresh,'Min',cfg.min)
set(handles.VolThresh,'Max',cfg.max)
set(handles.VolThresh,'SliderStep',[1/(cfg.max-cfg.min) 10/(cfg.max-cfg.min)])
set(handles.VolThresh,'Value',cfg.volThresh)
set(handles.VolThreshVal,'String',cfg.volThresh)

% register data
setappdata(handles.axesX,'mri',mri);
setappdata(handles.axesX,'cfg',cfg);
setappdata(handles.axesX,'subject',subject);
setappdata(handles.axesX,'focus',focus);
setappdata(handles.axesX,'iso',[]);

update_view(handles)
% UIWAIT makes mri_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mri_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderX_Callback(hObject, eventdata, handles)
% hObject    handle to sliderX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
focus=getappdata(handles.axesX,'focus');
focus(3)=round(get(hObject,'Value'));
setappdata(handles.axesX,'focus',focus);
update_view(handles)

% --- Executes during object creation, after setting all properties.
function sliderX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderY_Callback(hObject, eventdata, handles)
% hObject    handle to sliderY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
focus=getappdata(handles.axesX,'focus');
focus(2)=round(get(hObject,'Value'));
setappdata(handles.axesX,'focus',focus);
update_view(handles)

% --- Executes during object creation, after setting all properties.
function sliderY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
focus=getappdata(handles.axesX,'focus');
focus(1)=round(get(hObject,'Value'));
setappdata(handles.axesX,'focus',focus);
update_view(handles)

% --- Executes during object creation, after setting all properties.
function sliderZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function update_view(handles)
mri=getappdata(handles.axesX,'mri');
cfg=getappdata(handles.axesX,'cfg');
subject=getappdata(handles.axesX,'subject');
focus=getappdata(handles.axesX,'focus');

ax_items={'X','Y','Z'};

for d=1:length(ax_items)
 axes(handles.(['axes' ax_items{d}]));
 a=gca;
 switch d
  case 1
   slice=mri.anatomy(:,:,focus(3));
   vx=cfg.vx([1 2]);
   this_focus=focus([1 2]);
  case 2
   slice=reshape(mri.anatomy(:,focus(2),:),mri.dim(1),mri.dim(3));
   vx=cfg.vx([1 3]);
   this_focus=focus([1 3]);
   % z needs to be inverted
   this_focus(2)=mri.dim(3)-this_focus(2);
  case 3
   slice=reshape(mri.anatomy(focus(1),:,:),mri.dim(2),mri.dim(3));
   vx=cfg.vx([2 3]);
   this_focus=focus([2 3]);
   % z needs to be inverted
   this_focus(2)=mri.dim(3)-this_focus(2);
 end
 % deal with zoom
 if isfield(cfg,'zoom') && isnumeric(cfg.zoom) && cfg.zoom>0 
  full_fov=size(slice);
  zoomed_fov=round(full_fov/cfg.zoom);
  frame=[this_focus-round(zoomed_fov/2);this_focus+round(zoomed_fov/2)];
  if frame(1,1)<1
   frame(2,1)=frame(2,1)-frame(1,1);
   frame(1,1)=1;
  end
  if frame(1,2)<1
   frame(2,2)=frame(2,2)-frame(1,2);
   frame(1,2)=1;
  end
  if frame(2,1)>full_fov(1)
   frame(1,1)=frame(1,1)-(frame(2,1)-full_fov(1))+1;
   frame(2,1)=full_fov(1);
  end
  if frame(2,2)>full_fov(2)
   frame(1,2)=frame(1,2)-(frame(2,2)-full_fov(2))+1;
   frame(2,2)=full_fov(2);
  end
  slice=slice(frame(1,1):frame(2,1),frame(1,2):frame(2,2));
 end
 if vx(1)<0 
  slice=flipud(slice);
 end
 if vx(2)<0 
  slice=fliplr(slice);
 end
 
 % apply the windowing
 window_min=get(handles.WindowMin,'Value');
 window_max=get(handles.WindowMax,'Value');
 slice=((slice-window_min)*100)/(window_max-window_min);
 slice(slice>100)=100;
 slice(slice<0)=0;
 
 
 im=image(rot90(slice));
 a.YTickLabel=[];
 a.XTickLabel=[];
 im.ButtonDownFcn=@(hObject,eventdata)mri_view('FocusClick_Callback',hObject,eventdata,guidata(hObject),ax_items{d});
 % fix aspect ratio
 daspect([fliplr(abs(vx)) 1])
 colormap(gray)
 hold on
 % mark focus
 line([this_focus(1) this_focus(1)],[im.YData],'Color','w','LineWidth',0.5,'LineStyle','-.','HitTest','off','PickablePart','none');
 line([im.XData],[this_focus(2) this_focus(2)],'Color','w','LineWidth',0.5,'LineStyle','-.','HitTest','off','PickablePart','none');
 hold off
 
 % set the sliders
 set(handles.sliderX,'Value',focus(3))
 set(handles.sliderY,'Value',focus(2))
 set(handles.sliderZ,'Value',focus(1))

 %update 3D crosshair
 if isfield(cfg,'hz') && isvalid(cfg.hz)
  cfg.hz.YData=[focus(1) focus(1)];
  cfg.hz.XData=[focus(2) focus(2)];
  cfg.hz.ZData=[1 mri.dim(3)];
  
  cfg.hy.YData=[focus(1) focus(1)];
  cfg.hy.XData=[1 mri.dim(2)];
  cfg.hy.ZData=[focus(3) focus(3)]-mri.dim(3);

  cfg.hx.YData=[1 mri.dim(1)];
  cfg.hx.XData=[focus(2) focus(2)];
  cfg.hx.ZData=[focus(3) focus(3)]-mri.dim(3);

 end
end



%% helper functions - from SPM
function P = imatrix(M)
% returns the parameters for creating an affine transformation
% FORMAT P = spm_imatrix(M)
% M      - Affine transformation matrix
% P      - Parameters (see spm_matrix for definitions)
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Stefan Kiebel
% $Id$


% Translations and zooms
%-----------------------------------------------------------------------
R         = M(1:3,1:3);
C         = chol(R'*R);
P         = [M(1:3,4)' 0 0 0  diag(C)'  0 0 0];
if det(R)<0, P(7)=-P(7);end % Fix for -ve determinants

% Shears
%-----------------------------------------------------------------------
C         = diag(diag(C))\C;
P(10:12)  = C([4 7 8]);
R0        = my_matrix([0 0 0  0 0 0 P(7:12)]);
R0        = R0(1:3,1:3);
R1        = R/R0;

% This just leaves rotations in matrix R1
%-----------------------------------------------------------------------
%[          c5*c6,           c5*s6, s5]
%[-s4*s5*c6-c4*s6, -s4*s5*s6+c4*c6, s4*c5]
%[-c4*s5*c6+s4*s6, -c4*s5*s6-s4*c6, c4*c5]

P(5) = asin(rang(R1(1,3)));
if (abs(P(5))-pi/2)^2 < 1e-9,
    P(4) = 0;
    P(6) = atan2(-rang(R1(2,1)), rang(-R1(3,1)/R1(1,3)));
else
    c    = cos(P(5));
    P(4) = atan2(rang(R1(2,3)/c), rang(R1(3,3)/c));
    P(6) = atan2(rang(R1(1,2)/c), rang(R1(1,1)/c));
end;
return;

% There may be slight rounding errors making b>1 or b<-1.
function a = rang(b)
a = min(max(b, -1), 1);
return;

function [A] = my_matrix(P, order)
% returns an affine transformation matrix
% FORMAT [A] = spm_matrix(P, order)
% P(1)  - x translation
% P(2)  - y translation
% P(3)  - z translation
% P(4)  - x rotation about - {pitch} (radians)
% P(5)  - y rotation about - {roll}  (radians)
% P(6)  - z rotation about - {yaw}   (radians)
% P(7)  - x scaling
% P(8)  - y scaling
% P(9)  - z scaling
% P(10) - x affine
% P(11) - y affine
% P(12) - z affine
%
% order (optional) application order of transformations.
%
% A     - affine transformation matrix
%___________________________________________________________________________
%
% spm_matrix returns a matrix defining an orthogonal linear (translation,
% rotation, scaling or affine) transformation given a vector of
% parameters (P).  By default, the transformations are applied in the
% following order (i.e., the opposite to which they are specified):
%
% 1) shear
% 2) scale (zoom)
% 3) rotation - yaw, roll & pitch
% 4) translation
%
% This order can be changed by calling spm_matrix with a string as a
% second argument. This string may contain any valid MATLAB expression
% that returns a 4x4 matrix after evaluation. The special characters 'S',
% 'Z', 'R', 'T' can be used to reference the transformations 1)-4)
% above. The default order is 'T*R*Z*S', as described above.
%
% SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
% matrices of n coordinates.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id$


% pad P with 'null' parameters
%---------------------------------------------------------------------------
q  = [0 0 0 0 0 0 1 1 1 0 0 0];
P  = [P q((length(P) + 1):12)];

% default multiplication order if not specified
%---------------------------------------------------------------------------
if nargin < 2
    order = 'T*R*Z*S';
end;

T  =   [1   0   0   P(1);
        0   1   0   P(2);
        0   0   1   P(3);
        0   0   0   1];

R1  =  [1    0      0          0;
        0    cos(P(4))  sin(P(4))  0;
        0   -sin(P(4))  cos(P(4))  0;
        0    0      0          1];

R2  =  [cos(P(5))  0    sin(P(5))  0;
        0          1    0      0;
       -sin(P(5))  0    cos(P(5))  0;
        0          0    0          1];

R3  =  [cos(P(6))   sin(P(6))   0  0;
       -sin(P(6))   cos(P(6))   0  0;
        0           0           1  0;
        0           0       0  1];

R   = R1*R2*R3;

Z   =  [P(7)    0       0       0;
        0       P(8)    0       0;
        0       0       P(9)    0;
        0       0       0       1];

S   =  [1       P(10)   P(11)   0;
        0       1   P(12)   0;
        0       0       1   0;
        0       0       0       1];

A = eval(sprintf('%s;', order));
if ~isnumeric(A) || ndims(A) ~= 2 || any(size(A) ~= 4)
    error('Order expression ''%s'' did not return a valid 4x4 matrix.', ...
          order);
end;



function FocusClick_Callback(hObject, eventdata, handles, ax_item)
% hObject    handle to axesX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cfg=getappdata(handles.axesX,'cfg');
mri=getappdata(handles.axesX,'mri');
focus=getappdata(handles.axesX,'focus');

% now update the focus
switch ax_item
  case 'X'
   focus([1 2])=round(eventdata.IntersectionPoint([1 2]));

  case 'Y'
   focus([1 3])=round(eventdata.IntersectionPoint([1 2]));
   % z needs to be inverted
   focus(3)=mri.dim(3)-focus(3);
  case 'Z'
   focus([2 3])=round(eventdata.IntersectionPoint([1 2]));
   % z needs to be inverted
   focus(3)=mri.dim(3)-focus(3);
end
focus
setappdata(handles.axesX,'focus',focus);
%update
update_view(handles)


% --- Executes on slider movement.
function WindowMin_Callback(hObject, eventdata, handles)
% hObject    handle to WindowMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
window_min=get(handles.WindowMin,'Value');
set(handles.WindowMinVal,'String',window_min)
%update
update_view(handles)


% --- Executes during object creation, after setting all properties.
function WindowMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function WindowMax_Callback(hObject, eventdata, handles)
% hObject    handle to WindowMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

window_max=get(handles.WindowMax,'Value');
set(handles.WindowMaxVal,'String',window_max)
%update
update_view(handles)

% --- Executes during object creation, after setting all properties.
function WindowMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function WindowMinVal_Callback(hObject, eventdata, handles)
% hObject    handle to WindowMinVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowMinVal as text
%        str2double(get(hObject,'String')) returns contents of WindowMinVal as a double


% --- Executes during object creation, after setting all properties.
function WindowMinVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowMinVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WindowMaxVal_Callback(hObject, eventdata, handles)
% hObject    handle to WindowMaxVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowMaxVal as text
%        str2double(get(hObject,'String')) returns contents of WindowMaxVal as a double
volThresh=str2double(get(handles.VolThreshVal,'String'));
cfg=getappdata(handles.axesX,'cfg');
if volThresh>cfg.max
 volThresh=cfg.max;
end
if volThresh>cfg.min
 volThresh=cfg.min;
end
set(handles.VolThresh,'Value',volThresh)

% --- Executes during object creation, after setting all properties.
function WindowMaxVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowMaxVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function VolThresh_Callback(hObject, eventdata, handles)
% hObject    handle to VolThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.VolThreshVal,'String',get(handles.VolThresh,'Value'))


% --- Executes during object creation, after setting all properties.
function VolThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function VolThreshVal_Callback(hObject, eventdata, handles)
% hObject    handle to VolThreshVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VolThreshVal as text
%        str2double(get(hObject,'String')) returns contents of VolThreshVal as a double


% --- Executes during object creation, after setting all properties.
function VolThreshVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolThreshVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in VolEnable.
function VolEnable_Callback(hObject, eventdata, handles)
% hObject    handle to VolEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VolEnable

if handles.VolEnable.Value==1
 iso=getappdata(handles.axesX,'iso');
 % update 3D
 if isempty(iso)
  makeVol(handles);
 end
 updateVol(handles);
else
 % disable 3D
 cla(handles.axesVol);
 rotate3d('off')
 handles.VolRotate.Enable='off';
 handles.VolRotate.Value=0;
 handles.VolViewReset.Enable='off';
end


% update the 3D view
function updateVol(handles)
mri=getappdata(handles.axesX,'mri');
cfg=getappdata(handles.axesX,'cfg');
subject=getappdata(handles.axesX,'subject');
focus=getappdata(handles.axesX,'focus');

if handles.VolEnable.Value==1
 iso=getappdata(handles.axesX,'iso');
 if isempty(iso)
  makeVol(handles);
 end
 
 %get the right axes
 axes(handles.axesVol);
 a=gca;
 cfg.hiso=patch(a,iso,'Facecolor',[1,.75,.65],'EdgeColor','none');
 view(a,35,30)
 axis(a,'tight')
 axis off; % ticks and axes off
 axis(a,'vis3d'); % 3D rotatable
 lightangle(45,30)
 daspect(a,1./abs(cfg.vx))
 set(cfg.hiso, 'FaceLighting',   'gouraud');
 
 % draw focus lines
 hold on
 cfg.hz=plot3(a,[focus(1) focus(1)],[focus(2) focus(2)], [1 mri.dim(3)], 'Color','k','LineWidth',2,'HitTest','off','PickablePart','none');
 cfg.hy=plot3(a,[focus(1) focus(1)],[1 mri.dim(2)], [focus(3) focus(3)],  'Color','k','LineWidth',2,'HitTest','off','PickablePart','none');
 cfg.hx=plot3(a,[1 mri.dim(1)],[focus(2) focus(2)], [focus(3) focus(3)],  'Color','k','LineWidth',2,'HitTest','off','PickablePart','none');

 %uistack(cfg.hiso,'top');
 
 setappdata(handles.axesX,'cfg',cfg);
 handles.VolRotate.Enable='on';
 handles.VolViewReset.Enable='on';
end




% make a new 3d object
function makeVol(handles)
mri=getappdata(handles.axesX,'mri');
cfg=getappdata(handles.axesX,'cfg');
cfg.volThresh=handles.VolThresh.Value;

% make a new isosurface
msgbox('Generating a new 3D model, can take some time.')
iso=isosurface(smooth3(mri.anatomy),cfg.volThresh);
if length(iso.vertices)>500000
 iso=reducepatch(iso,500000);
end
handles.VolEnable.Enable='on';
setappdata(handles.axesX,'iso',iso);
setappdata(handles.axesX,'cfg',cfg);



function [norot]=myRotationFilter(obj,eventdata)
norot=false;
if isfield(get(obj),'ButtonDownFcn')
 norot= ~isempty(get(obj,'ButtonDownFcn'));
end


% --- Executes on button press in VolRotate.
function VolRotate_Callback(hObject, eventdata, handles)
% hObject    handle to VolRotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VolRotate
if handles.VolRotate.Value==1
 hrot=rotate3d;
 hrot.Enable='on';
 set(hrot,'ButtonDownFilter',@myRotationFilter);
else
 rotate3d('off')
end

% --- Executes on button press in VolRebuild.
function VolRebuild_Callback(hObject, eventdata, handles)
% hObject    handle to VolRebuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

makeVol(handles);
updateVol(handles);


% --- Executes on button press in VolViewReset.
function VolViewReset_Callback(hObject, eventdata, handles)
% hObject    handle to VolViewReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 axes(handles.axesVol);
 a=gca;
 view(a,35,30)
