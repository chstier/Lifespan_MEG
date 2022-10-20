function [ hFig ] = nmri_plot_surface_suma_interactive(surf, func, opt)
%[ hFig ] = nmri_plot_surface_suma_interactive(surf, func, opt)
%  
% This is the interavtive, network enabled variant of the SUMA surface
% viewer. It can deal with 2D network data and allows user selections of
% seeds
% 
% surf         =   surface structure,
%                  needs .pos, .tri and .curv
% func         =   functional map to overlay
% opt          =   options struct
%  .title      =   title
%  .output     =   output file
%  .scale      =   output scaling factor (increase for higher dpi)
%  .format     =   output file format, default PNG
%  .per_hemi   =   0 or 1 (split hemispheres)
%  .per_cortex =   0 or 1 (split cortex/subcortex)
%  .rot        =   rotation steps (azimuth, elevation pairs)
%              =   default [90 0; -90 0]
%  .rot_subcort=   rotation steps for subcortex (azimuth, elevation pairs)
%              =   default [90 60]
%  .zoom_subcort=  zoom faktor for subcortical, default 0.8
%  .roi_color  =   roi_lookup table
%  .thresh     =   binary abs threshold of func that is shown
%  .opathresh  =   opacity abs threshold of func that is shown
%  .fixalpha   =   opacity for all func vertices
%  .clim       =   colorscaling for funct [min max]
%  .colormap   =   colormap to use 

% written by NF 04/2017

% check the surface struct
if (~isfield(surf,'pos') || ~isfield(surf,'tri') || ~isfield(surf,'curv') )
 error('Could not find the necessary struct fields (.pos, .tri and .curv)')
end

% check func
if (~exist('func','var'))
 % make empty func
 func=zeros(length(surf.pos),1);
end

% make sure that func has the right orientation
if (size(func,1)>size(func,2))
 func=func';
end

% check length
if (length(surf.curv) ~= length(surf.pos))
 error('mismatch between curvature and vertices, should not happen')
end

% check dimensions
dim=length(size(func));
if dim>2
 error('Functional map needs to be 1D or 2D, more is not possible')
end

% check if 2D mode
if size(func,1)>1 && size(func,2)>1
 network_mode=true;
 netfunc=func;
 % make a condensed func to start with
 func=nanmean(func,1);
else
 netfunc=func;
 network_mode=false;
end

if (length(func) ~= length(surf.pos))
 error('mismatch between functional overlay and vertices, should not happen')
end

% check opt
if (~exist('opt','var'))
 % set defaults
 opt.per_hemi=0;
end

if (~isfield(opt,'per_hemi'))
 opt.per_hemi=0;
end

if (~isfield(opt,'per_cortex'))
 opt.per_cortex=0;
end

if (~isfield(opt,'rot'))
 if (opt.per_hemi>0)
  opt.rot = [-180 0; 0 0];
 else
  opt.rot = [180 0; -90 0; 0 0 ; 90 0];
 end
end

if (~isfield(opt,'rot_subcort'))
 opt.rot_subcort = opt.rot;
end

if (~isfield(opt,'zoom_subcort'))
 opt.zoom_subcort = 1;
end

if (~isfield(opt,'format')) 
 opt.format = 'png';
end

switch opt.format
 case 'png'
  print_format='-dpng';
  print_ext='png';
 case 'jpg'
  print_format='-djpeg';
  print_ext='jpg';
 case 'tiff'
  print_format='-dtiff';
  print_ext='tif';
 otherwise
  error('Print format not known')       
end

% check hemi split
if (isfield(opt,'per_hemi') && isfield(surf,'hemi') && opt.per_hemi==1)
 opt.nhemi=2;
else
 opt.nhemi=1;
end

if (isfield(opt,'per_cortex') && isfield(surf,'cortex') && opt.per_cortex==1)
 opt.ncortex=2;
else
 opt.ncortex=1;
end

% how many figuers will we get?
all_figs=opt.nhemi*opt.ncortex;
% set sizes
if all_figs>2 
 xsize=300*all_figs;
 ysize=400;
elseif all_figs==1
 % single
 xsize=800;
 ysize=600;
else
 % two
 xsize=900;
 ysize=500;
end 
if all_figs == 1
 ysize=xsize;
end
hFig=figure('Position',[0,0,xsize,ysize],'Visible','on');

% register data to Fig
setappdata(hFig,'netfunc',netfunc)
setappdata(hFig,'network_mode',network_mode)
setappdata(hFig,'opt',opt)
setappdata(hFig,'surf',surf)
setappdata(hFig,'func',func)
draw_brain(hFig)


end % end of function
 

function mouseCallback(obj,hit)
 hFig=get(get(obj,'Parent'),'Parent');
 netfunc=getappdata(hFig,'netfunc');
 surf=getappdata(hFig,'surf');
 network_mode=getappdata(hFig,'network_mode');
 % now determine closest vertex
 d=pdist2(double(hit.IntersectionPoint),double(obj.Vertices));
 [ds,iA]=sort(d);
 vertex=iA(1);
 
 % now find the correct vertex of all surf
 do=pdist2(double(obj.Vertices(vertex,:)), double(surf.pos));
 [dso,iAo]=sort(do);
 vertexo=iAo(1);
 
 % now check for network mode
 if network_mode
  % check is a legal vertex
  if ~isfield(surf,'annot') || surf.annot(vertexo)~=0
   % then redraw func on selection
   func=netfunc(vertexo,:);
   setappdata(hFig,'func',func);
   update_brain(hFig)
   col=get(obj,'FaceVertexCData');
   col(vertex,:)=repmat([0 1 0],length(vertex),1);
   set(obj,'FaceVertexCData',col);
  else
   % illegal vertex
   col=get(obj,'FaceVertexCData');
   col(vertex,:)=repmat([1 0 0],length(vertex),1);
   set(obj,'FaceVertexCData',col);
  end
 end
 if isfield(surf,'annot') 
  annot=nf_get_annotation(surf,vertexo);
  disp(sprintf('Vertex: %d, Annot: %s',vertexo,annot{1}))
 else
  disp(sprintf('Vertex: %d',vertexo))
 end
end


function draw_brain(hFig)
opt=getappdata(hFig,'opt');
surf=getappdata(hFig,'surf');
func=getappdata(hFig,'func');
% make color from curvature
cortex_light = [0.781 0.762 0.664];
cortex_dark  = [0.781 0.762 0.664]/2;
if (isfield(surf,'curv'))
 color = surf.curv(:) * cortex_dark + (1-surf.curv(:)) * cortex_light;
else
 color=[0.7 0.7 0.7];
end
% deal with ROI colors, if any
if (isfield(surf,'annot') && isfield(opt,'roi_color'))
 for i=1:length(opt.roi_color)
  color(surf.annot==opt.roi_color{i,1},:)=repmat(opt.roi_color{i,2},sum(surf.annot==opt.roi_color{i,1}),1);
 end
end

% deal with thresholds
if (isfield(opt,'thresh'))
 % hard binary threshold
 msk = double(abs(func)>opt.thresh)';
elseif (isfield(opt,'opathresh'))
 % soft threshold
 msk = double(abs(func)/opt.opathresh)';
 msk(msk>1)=1;
 msk(msk<0.1)=0;
elseif (isfield(opt,'fixalpha') && opt.fixalpha>=0 && opt.fixalpha<=1)
 msk=repmat(opt.fixalpha,size(func))';
else
 % make a somewhat meaningfull mask for all purposes
 msk = double(abs(func)/max(abs(func)))';
 msk(msk<0.3)=0;
end

% check annot for medial wall
if (isfield(surf,'annot'))
 % mask every func value without annotiation
 msk(surf.annot==0)=0;
 color(surf.annot==0,:)=repmat([0.6 0.6 0.6],sum(surf.annot==0),1);
end
% make inverse
msk_inv = 1-msk;
%clear figure, just in case
clf(hFig);
% deal with clim/caxis
if (isfield(opt,'clim') && length(opt.clim)==2)
 opt.cmin=opt.clim(1);
 opt.cmax=opt.clim(2);
elseif (sum(func<0)>0)
 opt.cmin=-max(func); 
 opt.cmax=max(func);
else
 opt.cmin=0; 
 opt.cmax=max(func);
end
fig_count=0;
all_figs=opt.nhemi*opt.ncortex;
for vcortex=1:opt.ncortex
 % loop for cortex
 for vhemi=1:opt.nhemi
  % loop per hemi
  fig_count=fig_count+1;
  % select the correct vertices
  sel_V=ones(length(surf.pos),1); % start with all
  sel_C=ones(length(surf.pos),1); % start with all

  % deal with alpha mask
  this_msk=msk;
  this_msk_inv=msk_inv;
  this_func=func;

  % check for cortex / subcortex
  add_brain=false;
  if (isfield(surf,'cortex') && opt.per_cortex==1)
   if (vcortex==1)
    sel_V(surf.cortex>1)=0;
   else
    % now we want to see subcortex
    this_msk(surf.cortex~=vcortex)=0; % set functional all to 0% except subcortex
    this_func(surf.cortex~=vcortex)=0; % set functional all to 0% except subcortex
    sel_C(surf.cortex>1)=0; % only keep cortex
    sel_V(surf.cortex<2)=0; % only keep subcortex
    add_brain=true;
   end
  end

  % check for hemi 
  if (isfield(surf,'hemi') && opt.per_hemi==1)
   sel_V(surf.hemi~=vhemi)=0;
   sel_C(surf.hemi~=vhemi)=0;
  end

  % make logical
  sel_V=sel_V>0;
  sel_C=sel_C>0;

  % now get hemi surface if needed
  this_surf=[];
  this_surf.pos=surf.pos;
  this_surf.tri=surf.tri;
  if (sum(sel_V)<length(surf.pos))
   % some need to go
   this_surf.pos(~sel_V,:)=repmat([NaN,NaN,NaN],sum(~sel_V),1);
   [this_surf.pos, this_surf.tri]=nf_remap_surface(this_surf.pos, this_surf.tri);
  end
  % generate a seperate brain - for subcortex
  if add_brain
   this_cortex=[];
   this_cortex.pos=surf.pos;
   this_cortex.tri=surf.tri;
   this_cortex.pos(~sel_C,:)=repmat([NaN,NaN,NaN],sum(~sel_C),1);
   [this_cortex.pos, this_cortex.tri]=nf_remap_surface(this_cortex.pos, this_cortex.tri);
  end

  if (all_figs>1)
   % we do not safe, so make subplots
   subplot(1,all_figs,fig_count)
   % and resize
   ax=gca;
   pos=ax.Position;
   pos(1)=0.01+((fig_count-1)/all_figs);
   pos(3)=(1/all_figs)*0.95;
   set(ax,'Position',pos)
  else
   clf(hFig)
  end

  this_col=color(sel_V,:); % brain color
  % set colormap 
  if (isfield(opt,'colormap'))
   switch opt.colormap
   case 'hot'
    colormap(hot(128))
   case 'cool'
    colormap(cool(128))
   case 'autumn'
    colormap(autumn(128))
   otherwise
    colormap('default') 
   end
  end
  cmap=colormap;

  % now map function values
  clength=size(cmap,1);
  % normalize func
  n_func=(this_func(sel_V)-opt.cmin)/(opt.cmax-opt.cmin);
  n_func(n_func>1)=1;
  n_func(n_func<0)=0;
  n_func=round(n_func*clength);
  % deal with NaN - set alpha to transparent
  this_msk(isnan(n_func))=0;
  this_msk_inv(isnan(n_func))=1;
  % and remove
  n_func(isnan(n_func))=1;
  % fix in range - just in case
  n_func(n_func<1)=1;
  n_func(n_func>clength)=clength;

  this_col=([color(sel_V,1).*this_msk_inv(sel_V) color(sel_V,2).*this_msk_inv(sel_V) color(sel_V,3).*this_msk_inv(sel_V)])+([cmap(n_func,1).*this_msk(sel_V) cmap(n_func,2).*this_msk(sel_V) cmap(n_func,3).*this_msk(sel_V)]);

  h1 =patch('Vertices', this_surf.pos, 'Faces', this_surf.tri, 'FaceVertexCData', this_col, 'FaceColor', 'interp','EdgeColor', 'none',...
      'ButtonDownFcn',@(obj,hit)mouseCallback(obj,hit));

  axis off; % ticks and axes off
  axis vis3d; % 3D rotatable
  axis equal; % scale to data
  caxis([opt.cmin opt.cmax]);
  set(h1, 'FaceLighting',   'gouraud');

  % add brain - for subcortex
  %if add_brain
  % h2=patch('Vertices', this_cortex.pos, 'Faces', this_cortex.tri, 'FaceVertexCData', color(sel_C,:), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.15,'HitTest','Off');
  %end
  
  % set start view
  if isfield(opt,'rot') && size(opt.rot,1)>1
   view(opt.rot(vhemi,:))
  else
   view([0 0])
  end
 end
end
end

function update_brain(hFig)
opt=getappdata(hFig,'opt');
surf=getappdata(hFig,'surf');
func=getappdata(hFig,'func');
% get the axes children
axs=get(hFig,'Children');
% make color from curvature
cortex_light = [0.781 0.762 0.664];
cortex_dark  = [0.781 0.762 0.664]/2;
if (isfield(surf,'curv'))
 color = surf.curv(:) * cortex_dark + (1-surf.curv(:)) * cortex_light;
else
 color=[0.7 0.7 0.7];
end
% deal with ROI colors, if any
if (isfield(surf,'annot') && isfield(opt,'roi_color'))
 for i=1:length(opt.roi_color)
  color(surf.annot==opt.roi_color{i,1},:)=repmat(opt.roi_color{i,2},sum(surf.annot==opt.roi_color{i,1}),1);
 end
end

% deal with thresholds
if (isfield(opt,'thresh'))
 % hard binary threshold
 msk = double(abs(func)>opt.thresh)';
elseif (isfield(opt,'opathresh'))
 % soft threshold
 msk = double(abs(func)/opt.opathresh)';
 msk(msk>1)=1;
 msk(msk<0.1)=0;
else
 % make a somewhat meaningfull mask for all purposes
 msk = double(abs(func)/max(abs(func)))';
 msk(msk<0.3)=0;
end
% check annot for medial wall
if (isfield(surf,'annot'))
 % mask every func value without annotiation
 msk(surf.annot==0)=0;
 color(surf.annot==0,:)=repmat([0.6 0.6 0.6],sum(surf.annot==0),1);
end
% make inverse
msk_inv = 1-msk;

% deal with clim/caxis
if (isfield(opt,'clim') && length(opt.clim)==2)
 opt.cmin=opt.clim(1);
 opt.cmax=opt.clim(2);
elseif (sum(func<0)>0)
 opt.cmin=-max(func); 
 opt.cmax=max(func);
else
 opt.cmin=0; 
 opt.cmax=max(func);
end
fig_count=0;
all_figs=opt.nhemi*opt.ncortex;
for vcortex=1:opt.ncortex
 % loop for cortex
 for vhemi=1:opt.nhemi
  % loop per hemi
  fig_count=fig_count+1;
  % select the correct vertices
  sel_V=ones(length(surf.pos),1); % start with all
  sel_C=ones(length(surf.pos),1); % start with all

  % deal with alpha mask
  this_msk=msk;
  this_msk_inv=msk_inv;
  this_func=func;

  % check for cortex / subcortex
  add_brain=false;
  if (isfield(surf,'cortex') && opt.per_cortex==1)
   if (vcortex==1)
    sel_V(surf.cortex>1)=0;
   else
    % now we want to see subcortex
    this_msk(surf.cortex~=vcortex)=0; % set functional all to 0% except subcortex
    this_func(surf.cortex~=vcortex)=0; % set functional all to 0% except subcortex
    sel_C(surf.cortex>1)=0; % only keep cortex
    sel_V(surf.cortex<2)=0; % only keep subcortex
    add_brain=true;
   end
  end

  % check for hemi 
  if (isfield(surf,'hemi') && opt.per_hemi==1)
   sel_V(surf.hemi~=vhemi)=0;
   sel_C(surf.hemi~=vhemi)=0;
  end

  % make logical
  sel_V=sel_V>0;
  sel_C=sel_C>0;

  % now get hemi surface if needed
  this_surf=[];
  this_surf.pos=surf.pos;
  this_surf.tri=surf.tri;
  if (sum(sel_V)<length(surf.pos))
   % some need to go
   this_surf.pos(~sel_V,:)=repmat([NaN,NaN,NaN],sum(~sel_V),1);
   [this_surf.pos, this_surf.tri]=nf_remap_surface(this_surf.pos, this_surf.tri);
  end
  % generate a seperate brain - for subcortex
  if add_brain
   this_cortex=[];
   this_cortex.pos=surf.pos;
   this_cortex.tri=surf.tri;
   this_cortex.pos(~sel_C,:)=repmat([NaN,NaN,NaN],sum(~sel_C),1);
   [this_cortex.pos, this_cortex.tri]=nf_remap_surface(this_cortex.pos, this_cortex.tri);
  end

  % make right axis active - somehow this is inverted order
  axes(axs(all_figs-fig_count+1));
  
  this_col=color(sel_V,:); % brain color
  % set colormap 
  if (isfield(opt,'colormap'))
   switch opt.colormap
   case 'hot'
    colormap(hot(128))
   case 'cool'
    colormap(cool(128))
   case 'autumn'
    colormap(autumn(128))
   otherwise
    colormap('default') 
   end
  end
  cmap=colormap;

  % now map function values
  clength=size(cmap,1);
  % normalize func
  n_func=(this_func(sel_V)-opt.cmin)/(opt.cmax-opt.cmin);
  n_func(n_func>1)=1;
  n_func(n_func<0)=0;
  n_func=round(n_func*clength);
  % deal with NaN - set alpha to transparent
  this_msk(sel_V(isnan(n_func)))=0;
  this_msk_inv(sel_V(isnan(n_func)))=1;
  % and remove
  n_func(isnan(n_func))=1;
  % fix in range - just in case
  n_func(n_func<1)=1;
  n_func(n_func>clength)=clength;

  this_col=([color(sel_V,1).*this_msk_inv(sel_V) color(sel_V,2).*this_msk_inv(sel_V) color(sel_V,3).*this_msk_inv(sel_V)])+([cmap(n_func,1).*this_msk(sel_V) cmap(n_func,2).*this_msk(sel_V) cmap(n_func,3).*this_msk(sel_V)]);

  % get h1
  patchh=get(axs(all_figs-fig_count+1),'Children');
  h1=patchh(end);
  set (h1,'FaceVertexCData',this_col);
  caxis([opt.cmin opt.cmax]);
 end
end



end