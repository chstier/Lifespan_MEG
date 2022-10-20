function [ choice ] = nf_uidialog( cfg )
%make our own dialog functions, maybe this can avoid MAtlab freezing the
%RDP interface on NMRI
% this one based on a figure NOT the dialog interface

    if ~exist('cfg','var')
     cfg=[];
    end
    if ~isfield(cfg,'title')
     cfg.title='Dialog Box';
    end
    if ~isfield(cfg,'question')
     cfg.question='Make your choice';
    end
    
    if ~isfield(cfg,'options')
     cfg.options={'Yes';'No'};
    end
    
    if ~isfield(cfg,'mode')
     cfg.mode='popup';
    end
    
    screensize = get(0, 'ScreenSize');
    win_pos=[0 0 300 150]; 
    if length(cfg.options)>2 && strcmpi('buttons',cfg.mode)
     % increase x for many buttons
     win_pos(3)=win_pos(3)+((length(cfg.options)-2)*80);
    end
    % find the length of the question
    txtl=max(cellfun(@strlen,cfg.question));
    if txtl>20 && win_pos(3)<(txtl*8)
     % increase x for long text
     win_pos(3)=max(win_pos(3),(txtl*8)+80);
    end
    
    if length(cfg.question)>1 
     % increase y for multi-line questions
     win_pos(4)=win_pos(4)+((length(cfg.question)-1)*18);
    end
    
    
    % center
    win_pos(1)=(screensize(3)-win_pos(3))/2;
    win_pos(2)=(screensize(4)-win_pos(4))/2;
    d=figure('Toolbar','none','DockControls','Off','MenuBar','none','Name',cfg.title,'NumberTitle','off',...
  'Position',win_pos,'Resize','off');
 
    
    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[(win_pos(3)/2)-((txtl*4)+5) 80 (txtl*8)+10 40+(length(cfg.question)*18)],...
           'String',cfg.question);
    
    if strcmpi('popup',cfg.mode)
     % have a popup choice list usefull for a lot of choices
     popl=max(cellfun(@strlen,cfg.options));
     popup = uicontrol('Parent',d,...
           'Style','popup',...
           'Position',[(win_pos(3)/2)-((popl*4)+5) 70 (popl*8)+10 25],...    
           'String',cfg.options,...
           'Callback',@popup_callback);
       
     confb = uicontrol('Parent',d,...
           'Position',[(win_pos(3)/2)-130 20 120 25],...
           'String','Select and Close',...
           'Callback', 'delete(gcf)');
          
     abortb = uicontrol('Parent',d,...
           'Position',[(win_pos(3)/2)+20 20 120 25],...
           'String','Abort',...
           'Callback',@abort_callback);
          
    elseif strcmpi('buttons',cfg.mode)
     % give choice buttons
     for i=1:length(cfg.options)
      xw=round(win_pos(3)/(length(cfg.options)+1));
      % check if default
      if isfield(cfg,'default') && strcmpi(cfg.default,cfg.options{i})
       fw='bold';
      else
       fw='normal';
      end
      
      btn{i} = uicontrol('Parent',d,...
           'Position',[(xw*(i-1))+20 20 120 25],...
           'String',cfg.options{i},...
           'FontWeight',fw,...
           'Callback',@btn_callback);
     end
     % add abort
       abortb = uicontrol('Parent',d,...
           'Position',[(xw*(i))+20 20 120 25],....
           'String','Abort',...
           'Callback',@abort_callback);
     
    end
    
    
    % default choice
    if isfield(cfg,'default')
     choice=cfg.default;
    else
     choice =  cfg.options{1};
    end
    
    % Wait for d to close before running to completion
    uiwait(d);
   
       function popup_callback(popup,event)
          idx = popup.Value;
          popup_items = popup.String;
          % This code uses dot notation to get properties.
          % Dot notation runs in R2014b and later.
          % For R2014a and earlier:
          % idx = get(popup,'Value');
          % popup_items = get(popup,'String');
          choice = char(popup_items(idx,:));
       end
      
       function abort_callback(popup,event)
          choice = '';
          delete(gcf)
       end
      
       function btn_callback(popup,event)
          choice = popup.String;
          delete(gcf)
       end
end


