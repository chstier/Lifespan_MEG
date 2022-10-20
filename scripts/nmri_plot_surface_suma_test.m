function [ hFig ] = nmri_plot_surface_suma_test(surf, func, opt)
%[ hFig ] = nmri_plot_surface(surf, func, opt)
%  
% Will make standard plots of SUMA derived surfaces 
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

% written by NF 11/2016

% check the surface struct
if (~isfield(surf,'pos') || ~isfield(surf,'tri') || ~isfield(surf,'curv') )
 error('Could not find the necessary struct fields (.pos, .tri and .curv)')
end

% check func
if (~exist('func','var')) || isempty(func)
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
 
% check length
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
  opt.rot = [0 0; 180 0];
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


% check hemi split
if (isfield(opt,'per_hemi') && isfield(surf,'hemi') && opt.per_hemi==1)
 nhemi=2;
else
 nhemi=1;
end

if (isfield(opt,'per_cortex') && isfield(surf,'cortex') && opt.per_cortex==1)
 ncortex=2;
else
 ncortex=1;
end

% how many figuers will we get?
if isfield(opt,'output')
 all_figs=size(opt.rot,1)*nhemi;
 if ncortex==2
  all_figs=all_figs+(size(opt.rot_subcort,1)*nhemi);
 end
else
 all_figs=nhemi*ncortex;
end


% set sizes
if (~isfield(opt,'output'))
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
 if all_figs > 1

 else
  ysize=xsize;
 end
 hFig=figure('Position',[0,0,xsize,ysize],'Visible','on');
else
 imgX=800;
 imgY=800;
 if (~isfield(opt,'scale'))
  opt.scale=1;
 end
 print_dpi=['-r' num2str(200*opt.scale)];
 print_dpi_col=['-r' num2str(150*opt.scale)];
 hFig=figure('Position',[0,0,imgX,imgY],'Visible','off');
end
fig_count=1;
all_plots='';
plot_list={};
plot_cortex={};


    
%clear figure, just in case
clf(hFig);

% deal with output info
if (isfield(opt,'output'))
 [out_pa out_file out_ext]=fileparts(opt.output);
 if (~exist(out_pa,'dir') && ~isempty(out_pa))
  mkdir(out_pa)
 end
 out_ext=['.' print_ext];
end

% set export for printer
set(hFig,'PaperUnits','inches','PaperPosition',[0 0 4 4]) % set square output


% deal with clim/caxis
if (isfield(opt,'clim') && length(opt.clim)==2)
 cmin=opt.clim(1);
 cmax=opt.clim(2);
elseif (sum(func<0)>0)
 cmin=-max(func); 
 cmax=max(func);
else
 cmin=0; 
 cmax=max(func);
end

for vcortex=1:ncortex
 % loop for cortex
 for vhemi=1:nhemi
  % loop per hemi
   
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
  
  
  if (~isfield(opt,'output') && (nhemi*ncortex>1))
   % we do not safe, so make subplots
   subplot(1,nhemi*ncortex,fig_count)
  else
   clf(hFig)
  end
 
  this_col=color(sel_V,:); % brain color
  % set colormap 
  if (isfield(opt,'colormap'))
   switch opt.colormap
   case 'hot'
    colormap(hot(256))
   case 'cool'
    colormap(cool(128))
   case 'autumn'
    colormap(autumn(128))
   case 'jet'
    colormap(jet(128))
   otherwise
    colormap(parula(128)) 
   end
  end
  cmap=colormap;
  
  % now map function values
  clength=size(cmap,1);
  % normalize func
  n_func=(this_func(sel_V)-cmin)/(cmax-cmin);
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
    
  this_col=([color(sel_V,1).*this_msk_inv(sel_V) color(sel_V,2).*this_msk_inv(sel_V) color(sel_V,3).*this_msk_inv(sel_V)])...
   +([cmap(n_func,1).*this_msk(sel_V) cmap(n_func,2).*this_msk(sel_V) cmap(n_func,3).*this_msk(sel_V)]);
  
  h1 =patch('Vertices', this_surf.pos, 'Faces', this_surf.tri, 'FaceVertexCData', this_col, 'FaceColor', 'interp','EdgeColor', 'none');
  axis off; % ticks and axes off
  axis vis3d; % 3D rotatable
  axis equal; % scale to data
  caxis([cmin cmax]);

  set(h1, 'FaceLighting',   'gouraud');

  % add brain - for subcortex
  if add_brain
   h2=patch('Vertices', this_cortex.pos, 'Faces', this_cortex.tri, 'FaceVertexCData', color(sel_C,:), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
  end
  
  if (isfield(opt,'output'))
   % save images
   % check rotation
   if (vcortex==2)
    this_rot=opt.rot_subcort;
   else
    this_rot=opt.rot;
   end
   for Ri=1:size(this_rot,1)
    view(this_rot(Ri,:))
    this_plot=fullfile(out_pa,[out_file '_suma_subplot_' sprintf('%02d',fig_count) '.' print_ext]);
    all_plots=[all_plots ' ' this_plot];
    plot_list{end+1}=this_plot;
    plot_cortex{end+1}=vcortex;
    print(hFig,this_plot,print_format,print_dpi)
    if fig_count==all_figs
     % add colorbar and save as well   
     view([0 0])
     % hide patches
     if exist('h1','var')
      set(h1,'Visible','Off')
     end
     if exist('h2','var')
      set(h2,'Visible','Off')
     end
     colorbar
     colbar=fullfile(out_pa,[out_file '_suma_subplot_col.' print_ext]);
     print(hFig,colbar,print_format,print_dpi_col)
     % only show after last save, if just one fig and not saved
     if ~isfield(opt,'output') || all_figs==1 
      if usejava('desktop') 
       % do not make visible if in non-java mode
       set(hFig,'Visible','on')
       if exist('h1','var')
        set(h1,'Visible','On')
       end
       if exist('h2','var')
        set(h2,'Visible','On')
       end
      else
       warning('Would not make plot visible in non-JAVA enabled mode')
      end
     end
    else
     fig_count=fig_count+1;
    end
   end
  else
   % no save, so show
   % make sensible roation and zoom
   if (isfield(surf,'hemi') && opt.per_hemi==1)
    if vhemi==1
     view([-90 0]);
    else
     view([90 0]);
    end
   else
    view([180 0]);
   end
   if fig_count==(nhemi*ncortex)
    % add colorbar, if single plot
    if all_figs==1
     colorbar('EastOutside')
    end
    if usejava('desktop')
     % do not make visible if in non-java mode
     set(hFig,'Visible','on')
    else
     warning('Would not make plot visible in non-JAVA enabled mode')
    end
   else
    fig_count=fig_count+1;
   end
  end
 end
end


if (isfield(opt,'output'))
 % now do the montage of the plots

 if (exist('colbar','var'))
  % start by cutting the colorbar
  colbar_clip=fullfile(out_pa,[out_file '_suma_subplot_col_clip.' print_ext]);
  cmd=['convert ' colbar ' -crop ' sprintf('%dx%d+%d+%d',180*opt.scale,imgY*opt.scale,420*opt.scale,0) ' ' colbar_clip];
  system(cmd);
 else
  error('no colorbar image found--should not happen')
  colbar_clip='';
 end
 
 % apply some y clipping to images and makes halves
 plots_left=[];
 plots_right=[];
 
 for n=1:length(plot_list)
  cmd=['convert ' plot_list{n} ' -crop ' sprintf('%dx%d+%d+%d',(imgX*opt.scale)*0.9,(imgY*opt.scale)*0.8,(imgX*opt.scale)*0.05,(imgY*opt.scale)*0.1) ' ' plot_list{n}];
  system(cmd);
  if n<=(round(length(plot_list)/2))
   if isempty(plots_left)
     plots_left=plot_list{n};
   else
    plots_left=[plots_left ' ' plot_list{n}];
   end
  else
   if isempty(plots_right)
    plots_right=plot_list{n};
   else
    plots_right=[plots_right ' ' plot_list{n}];
   end
  end
 end
 
 output_all=fullfile(out_pa,[out_file out_ext]);
 % merge all plots
 cmd=['montage -gravity South -mode concatenate -tile ' num2str(fig_count+1) 'x -background white ' plots_left ' ' colbar_clip ' ' plots_right ' ' output_all];
 system(cmd);
 
 % and annotate if wanted
 if (isfield(opt,'title'))
  cmd=['convert ' output_all ' -font /usr/share/fonts/liberation/LiberationSans-Bold.ttf -pointsize ' num2str(36*opt.scale) ' -fill black -background white -gravity North -splice 0x48 -annotate +0+2 "' opt.title '" ' output_all];
  system(cmd);
 end
 
 % delete temporary
 system(['rm -f ' all_plots ' ' colbar ' ' colbar_clip]);
end

end % end of function
 
%view([90 0]) % medial view (lh) - lateral view (rh)
% view([180 0]) % ap view
% view([-90 0]) % lateral view (lh) - medial view (rh)
% view([0 0]) % post view

