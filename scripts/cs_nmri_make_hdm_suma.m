function [ subject, data ] = cs_nmri_make_hdm_suma(subject, data, params)
%[ subject, data ] = nmri_make_hdm_suma(subject, data, params)
%  
% This function will deal with reading the MRI info 
% per default will read from SUMA, other options we be integrated later 
% will assume that the subject.id = basename
% otherwise you should specify subject.suma_dir
% 
% subject   =   subject structure, see nmri_read_subject
%               could be a struct or .m file or .mat file
% data      =   will return the data (optional)

% written by NF 11/2016 and modified by Christina Stier, 2020/22


% Freedurder storage dir
fsdir = getenv('SUBJECTS_DIR'); % if Freesurfer is correctly setup, this shoudl point to the right place


% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

% check if we have clean_dataset (maybe ICA or not)
if (isfield(params,'useICA_clean') && params.useICA_clean==1)
  if (isfield(subject,'cleanICA_dataset') &&  exist(subject.cleanICA_dataset,'file'))
   input=subject.cleanICA_dataset;
  else
   error('ICA cleaned dataset not found - Run nmri_artifactrejection first')
  end
else
%  if (isfield(subject,'clean_dataset') &&  exist(subject.clean_dataset,'file'))
%    input=subject.clean_dataset;
%   else
%    error('Cleaned dataset (w/o ICA) not found - Run nmri_artifactrejection first')
%  end  
end

if (~exist('data','var') || isempty(data) )  
 disp(sprintf('Loading data = %s ...',input))
 load(input,'data');
 
 if (~exist('data','var') ) 
  error('Could not load data')
 end
end

% call central dataset selector
subject=nmri_determine_datatype(subject);

%% load dws-dataset to get grad-info
prev_subj = load(subject.dws_filt_dataset);
data.grad = prev_subj.data.grad;

%% Get the modality-specific analysis params
[ params ] = nmri_get_modality_params( params, subject.dtype );


%% make our output paths and dir

if (~isfield(subject,'mri_dir'))
 subject.mri_dir=fullfile(subject.analysis_dir,subject.id,'mri');
end
if (~exist(subject.mri_dir,'dir'))
 mkdir(subject.mri_dir)
end

if (~isfield(subject,'proc_dir'))
 subject.proc_dir=fullfile(subject.analysis_dir,subject.id,'processed');
end
if (~exist(subject.proc_dir,'dir'))
 mkdir(subject.proc_dir)
end

% These are subject specific, so could be re-used for all exams
if (~isfield(subject,'mri_neuromag'))
 subject.mri_neuromag=fullfile(subject.mri_dir,['mri_neuromag_' subject.id '.mat']);
end
if (~isfield(subject,'mri_seg'))
 subject.mri_seg=fullfile(subject.mri_dir,['mri_seg_' subject.id '.mat']);
end
if (~isfield(subject,'mri_T1'))
 subject.mri_T1=fullfile(subject.mri_dir,['mri_T1_' subject.id '.nii']);
end
if (~isfield(subject,'suma_surface'))
 subject.suma_surface=fullfile(subject.mri_dir,['suma_surface_' subject.id '.mat']);
end

% This is subject and exam_id specific (leadfield may differ)
if (~isfield(subject,'hdm_lead'))
 % always set to current standard
 subject.hdm_lead=fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_' params.headmodel '.mat']);
end

%% check the SUMA and segmentation and run if needed
if (~isfield(subject,'suma_dir') && ( ~exist(subject.mri_neuromag,'file') || ~exist(subject.mri_seg,'file') || ~exist(subject.suma_surface,'file')))
 % not have the needed info - call the mri suma script to deal with that
 subject=cs_nmri_mri_suma_segment(subject);
end

if (~exist(subject.mri_neuromag,'file') || ~exist(subject.mri_seg,'file') || ~exist(subject.suma_surface,'file'))
 % not have the needed info - call the mri suma script to deal with that
 subject=cs_nmri_mri_suma_segment(subject);
end

%% Read MRI - interactive part

% Load T1 and as to manually mark LPA - RPA - nasion and Z-dir
if (exist(subject.mri_neuromag,'file'))
 disp('Found neuromag-aligend MRI - loading')
 load(subject.mri_neuromag,'mri_neuromag','mri')
else
 % not present - should not happen
 error('Problems with getting the neuromag-aligned MRI')
end
  
%% Load MRI seg

% check the MRI segmentation
if (exist(subject.mri_seg,'file'))
 disp('Found segmented MRI - loading')
 load(subject.mri_seg,'mri_seg')
else
 % not present - should not happen
 error('Problems with getting the segmented MRI')
end

%% now comes the SUMA / ASEG mapping
if (exist(subject.suma_surface,'file'))
 disp('Found SUMA surface - loading')
 load(subject.suma_surface,'suma_all');
else
 % not present - should not happen
 error('Problems with getting the SUMA surface')
end

% Transform mri_seg to neuromag, if needed
if (~isfield(mri_seg,'coordsys') || ~strcmpi(mri_seg.coordsys,'neuromag'))
 % note the transformation difference of manual alignement
 T = mri_neuromag.transform*inv(mri_neuromag.transformorig);
 
 mri_seg = ft_convert_units(mri_seg, 'cm');
 mri_seg = ft_transform_geometry(T, mri_seg);
 mri_seg.coordsys='neuromag';
end


%% Meshing Brain - needed for singleshell models and visualization
cfg           = [];
cfg.tissue    = {'brain','skull','scalp'};
cfg.method    = 'projectmesh';
cfg.numvertices = [3000 2000 3000]; 

% make mri_bnd for brain and scalp 
mri_bnd = [];
mri_bnd.transform = mri_seg.transform;
mri_bnd.brain = mri_seg.brain;
% make sure to have all SUMA points as brain
%for i=1:length(suma_all.pos)
% xyz=round(inv(mri_seg.transform)*[suma_all.pos(i,:) 1]');
% mri_bnd.brain(xyz(1:3))=1;
%end
braindil = imdilate(mri_seg.brain>0, strel_bol(6)); % from fieldtrip ft_volumesegment
mri_bnd.skull = braindil > 0.5 ;
mri_bnd.scalp = mri_seg.tissue;
mri_bnd.dim = mri_seg.dim;
mri_bnd.coordsys = mri_seg.coordsys;
mri_bnd.unit = mri_seg.unit;

bnd = ft_prepare_mesh(cfg,mri_bnd);

clear mri_bnd braindil

% smooth the brain (bnd(1)) mesh with iso2mesh function
%bnd(1).pos=smoothsurf(bnd(1).pos,[],meshconn(bnd(1).tri,length(bnd(1).pos)),1,1,'lowpass');
%bnd(2).pos=smoothsurf(bnd(2).pos,[],meshconn(bnd(2).tri,length(bnd(2).pos)),1,1,'lowpass');
% smooth the scalp (bnd(3)) mesh with iso2mesh function - more aggressive
bnd(3).pos=smoothsurf(bnd(3).pos,[],meshconn(bnd(3).tri,length(bnd(3).pos)),1,1,'lowpass');


%% Now do the electrode alignement for EEG
if (strcmp(subject.dtype,'EEG'))
 if (~isfield(subject,'electrodes_aligned'))
  subject.electrodes_aligned=fullfile(subject.analysis_dir,subject.id,'processed',['electrodes_aligned_' subject.id '_' subject.exam_id '.mat']);
 end
 
 if (exist(subject.electrodes_aligned,'file'))
  disp('Found aligend EEG electrodes - loading')
  load(subject.electrodes_aligned,'elec_aligned','elec_present','elec_missing')
 else
  if (~isfield(subject,'elec_file') || ~isfield(subject,'elec_fiducials') )
   error('Electrodes information missing in suject struct, should be auto-determined by dataset selector. Check dataset-mappings files.')
  end
  
  % read and re-align sensors
  sens = ft_read_sens(subject.elec_file);
  
  geoscan_txt='';
  
  % check for Geoscan (EGI 3D electrode pos scan)
  if exist(fullfile(subject.raw_dataset,'geoscan.txt'),'file') 
   if params.Geoscan_Use 
    % update by the Geoscan
    fid=fopen(fullfile(subject.raw_dataset,'geoscan.txt'),'r');
    out=textscan(fid,'%s %f %f %f %d','CommentStyle','#','Headerlines',5);
    fclose(fid);
    
    % now loop for the channels
    for i=1:length(sens.label)
     ind=find(strcmp(sens.label{i},out{1}));
     if ~isempty(ind)
      % found this one - update
      sens.chanpos(i,:)=[out{2}(ind) out{3}(ind) out{4}(ind)];  
      sens.elecpos(i,:)=[out{2}(ind) out{3}(ind) out{4}(ind)];
     else
      error('Could not find electrode=%s in Geoscan file=%s. Investigate and re-generate the file if needed.',sens.label{i},fullfile(subject.raw_dataset,'geoscan.txt'))
     end
    end
   end
   geoscan_txt='[using Geoscan data]';
  elseif params.Geoscan_Require
   error('A Geoscan is required (as per analysis_params), but it was not found in the raw dataset (expected=%s)',fullfile(subject.raw_dataset,'geoscan.txt'))
  end
  
  
            
  % fiducials from MRI
  nas = mri_neuromag.cfg.fiducial.nas;
  lpa = mri_neuromag.cfg.fiducial.lpa;
  rpa = mri_neuromag.cfg.fiducial.rpa;
  nas=ft_warp_apply(mri_neuromag.transform,nas, 'homogenous');
  lpa=ft_warp_apply(mri_neuromag.transform,lpa, 'homogenous');
  rpa=ft_warp_apply(mri_neuromag.transform,rpa, 'homogenous');
        
  % create a structure similar to a template set of electrodes - updated
  % for target
  fid=[];
  fid.elecpos       = [nas; lpa; rpa];   
  fid.chanpos       = [nas; lpa; rpa]; % neuromag-coordinates of fiducials               
  fid.label         = subject.elec_fiducials;    % same labels as in elec
         
  % alignment - First by fiducial points
  cfg               = [];
  cfg.method        = 'fiducial';
  cfg.target        = fid;                   % see above
  cfg.elec          = sens;
  
  cfg.fiducial      = subject.elec_fiducials;  % labels of fiducials in fid and in elec
  elec_aligned      = ft_electroderealign(cfg);

  if (isfield(params,'QC_plots') && params.QC_plots==1)
   % make the QC dir
   QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
   if (~exist(QCdir,'dir'))
    mkdir(QCdir)
   end
 
   hFig=figure('Position',[0,0,400,400],'Visible','off');
   hold on;
   ft_plot_sens(elec_aligned,'style','.','elecsize',10);
   ft_plot_mesh(bnd(3),'facecolor',[1 0.8 0.6],'edgecolor','none','facealpha',0.8);
   ft_plot_mesh(nas(1:3),'vertexcolor','w','vertexsize',20,'vertexmarker','x');
   ft_plot_mesh(lpa(1:3),'vertexcolor',[0 1 0],'vertexsize',20,'vertexmarker','x');
   ft_plot_mesh(rpa(1:3),'vertexcolor',[1 0 0],'vertexsize',20,'vertexmarker','x');
   camlight('headlight'); 
   hold off;
   hFig=nmri_rotate_print(hFig,[0 0;180 0;90 0;0 90],fullfile(QCdir,['electrodes_fiducials_' subject.id '_' subject.exam_id]),[subject.id ' - Fiducials only' geoscan_txt]);
   % only show after save
   title([subject.id ' - Fiducials only' geoscan_txt]);
   view([0 0])
   set(hFig,'Visible','on')
  end
  
  % now project on headshape
  cfg           = [];
  cfg.method    = 'project';
  cfg.elec      = elec_aligned;
  cfg.headshape = bnd(3);
  elec_aligned  = ft_electroderealign(cfg);
  elec_aligned.label{end} = 'VREF';
  
  % and plot again
  if (isfield(params,'QC_plots') && params.QC_plots==1)
   hFig=figure('Position',[0,0,400,400],'Visible','off');
   hold on;
   ft_plot_sens(elec_aligned,'style','.','elecsize',10);
   ft_plot_mesh(bnd(3),'facecolor','skin','edgecolor','none','facealpha',0.8);
   ft_plot_mesh(nas(1:3),'vertexcolor','w','vertexsize',20,'vertexmarker','x');
   ft_plot_mesh(lpa(1:3),'vertexcolor',[0 1 0],'vertexsize',20,'vertexmarker','x');
   ft_plot_mesh(rpa(1:3),'vertexcolor',[1 0 0],'vertexsize',20,'vertexmarker','x');
   camlight('headlight');
   hold off;

   hFig=nmri_rotate_print(hFig,[0 0;180 0;90 0;0 90],fullfile(QCdir,['electrodes_projected_' subject.id '_' subject.exam_id]),[subject.id ' - Projected' geoscan_txt]);
   % only show after save
   title([subject.id ' - Projected' geoscan_txt]);
   view([0 0])
   set(hFig,'Visible','on')
  end 
 
  % sort electrodes into missing and present
  elec_present=[];
  elec_present.unit=elec_aligned.unit;
  elec_present.label={};
  elec_present.chantype={};
  elec_present.chanunit={};
  elec_present.chanpos=[];
  elec_present.elecpos=[];
  elec_missing=[];
  elec_missing.unit=elec_aligned.unit;
  elec_missing.label={};
  elec_missing.chantype={};
  elec_missing.chanunit={};
  elec_missing.chanpos=[]; 
  elec_missing.elecpos=[]; 
  % make UPPER case labels -- seems important
  elec_aligned.label=upper(elec_aligned.label);
  for i=1:length(elec_aligned.label)
   match_elec=find(strcmp(data.label,elec_aligned.label{i})); % NOT do case insens.
   
   % deal with channels marked as BAD, as of 10/2019
   if isfield(data,'bad_channels') && any(strcmp(elec_aligned.label{i},data.bad_channels))
    % electrode is present, but bad  -treat as missing
    match_elec=[];
   end
   
   if isempty(match_elec)
    elec_missing.label{end+1,1}=elec_aligned.label{i};
    elec_missing.chantype{end+1,1}=elec_aligned.chantype{i,1};
    elec_missing.chanunit{end+1,1}=elec_aligned.chanunit{i,1};
    elec_missing.chanpos(end+1,1:3)=elec_aligned.chanpos(i,1:3); 
    elec_missing.elecpos(end+1,1:3)=elec_aligned.elecpos(i,1:3); 
   else
    elec_present.label{end+1,1}=elec_aligned.label{i,1};
    elec_present.chantype{end+1,1}=elec_aligned.chantype{i,1};
    elec_present.chanunit{end+1,1}=elec_aligned.chanunit{i,1};
    elec_present.chanpos(end+1,1:3)=elec_aligned.chanpos(i,1:3); 
    elec_present.elecpos(end+1,1:3)=elec_aligned.elecpos(i,1:3);
   end
  end
   
  save(subject.electrodes_aligned,'elec_aligned','elec_present','elec_missing');
 end
end

%% Check if we have a headmodel already, maybe from a different exam_id
check_hdm=dir(fullfile(subject.proc_dir,['hdm_lead_' subject.id '_*_' params.headmodel '.mat']));
if size(check_hdm,1)==1
 % we have found a single headmodel, so load
 fprintf('Found a previously calculated headmodel (%s) - re-using\n\n',check_hdm(1).name)
 load(fullfile(subject.proc_dir,check_hdm(1).name),'hdm','bnd')
end 

 
%% Build FEM model, if needed
if (strcmp(params.headmodel,'simbio') && ~exist('hdm','var'))
 % simbio does not like Freesurfe-style transformations - hence we need to 
 % resample all classes to more standard transform
 cfg=[];
 cfg.dim=mri_seg.dim;
 mri_bnd = [];
 mri_bnd.transform = mri_seg.transform;
 mri_bnd.dim = mri_seg.dim;
 mri_bnd.coordsys = mri_seg.coordsys;
 mri_bnd.unit = mri_seg.unit;
 all_fields={'grey','white','csf','skull','tissue'};
 for i=1:length(all_fields)
  mri_this=mri_bnd;
  mri_this.anatomy=double(mri_seg.(all_fields{i}));
  mri_this=ft_convert_units(mri_this,'mm');
  mri_this=ft_volumereslice(cfg,mri_this);
  mri_this=ft_convert_units(mri_this,'cm');
  mri_bnd.(all_fields{i})=mri_this.anatomy>0.5; % rebinarize
 end
 
 % we want to clean the classes with the smooth scalp outline from above
 scalp_vox=ft_transform_geometry(inv(mri_seg.transform),bnd(3));
 FV.vertices=scalp_vox.pos;
 FV.faces=scalp_vox.tri;
 scalp_img=polygon2voxel(FV,mri_seg.dim,'none',0);
 a1 = nf_volumefillholes_all(scalp_img, 1);
 a2 = nf_volumefillholes_all(scalp_img, 2);
 a3 = nf_volumefillholes_all(scalp_img, 3);
 scalp_img=a1|a2|a3;
 
 % now also transfrom scalp
 mri_this=mri_bnd;
 mri_this.anatomy=double(scalp_img);
 mri_this=ft_convert_units(mri_this,'mm');
 mri_this=ft_volumereslice(cfg,mri_this);
 mri_this=ft_convert_units(mri_this,'cm');
 scalp_img=mri_this.anatomy>0.5; % rebinarize
  
 
 % now update mri_bnd
 mri_bnd.transform = mri_this.transform;
 mri_bnd.dim = mri_this.dim;
 mri_bnd.coordsys = mri_this.coordsys;
 mri_bnd.unit = mri_this.unit;
 % deal with tissue/scalp
 mri_bnd.scalp=mri_bnd.tissue;
 mri_bnd=rmfield(mri_bnd,'tissue');
  
 % now mask
 mri_bnd.grey(scalp_img==0)=0;
 mri_bnd.white(scalp_img==0)=0;
 mri_bnd.csf(scalp_img==0)=0;
 mri_bnd.skull(scalp_img==0)=0;
 mri_bnd.scalp(scalp_img==0)=0; 
  
 cfg        = [];
 cfg.shift  = 0.3;
 cfg.method = 'hexahedral';
 fem = ft_prepare_mesh(cfg,mri_bnd);
 clear mri_bnd scalp_img a1 a2 a3
end

%% Headmodeling
if (~exist('hdm','var'))
 cfg          = [];
 thdm         = [];
 if (strcmp(subject.dtype,'MEG'))
  % MEG  
  switch params.headmodel
  case 'singleshell'
   cfg.method    = 'singleshell';
   if (isfield(params,'headmodel_class') && strcmp(params.headmodel_class,'suma'))
    % use SUMA mesh directly
    thdm.pos     = suma_all.pos;
    thdm.tri     = suma_all.tri;
   else
   % default to use brain (from segemented mri)
    thdm.pos     = bnd(1).pos;
    thdm.tri     = bnd(1).tri;
   end
   otherwise
   error('unsupported head model for MEG')
  end
 elseif (strcmp(subject.dtype,'EEG'))
  % EEG
  switch params.headmodel
  case 'dipoli'
   % BEM head model
   cfg.method    = 'dipoli';
   % downsample scalp / bnd(3) to 1000 , standard in fieldtrip and faster,
   target_vertices=1000;
   target_faces=target_vertices*2; % some first target
   this_mesh.faces=bnd(3).tri;
   this_mesh.vertices=bnd(3).pos;
   while(length(this_mesh.vertices)>target_vertices)
    [this_mesh.faces, this_mesh.vertices]=reducepatch(this_mesh,target_faces);
    % reduce by further 2% if needed
    target_faces=floor(length(this_mesh.faces)*0.98);
    fprintf('Reducing scalp mesh (vx count=%d)\n',length(this_mesh.vertices))
   end
   bnd(3).pos=this_mesh.vertices;
   bnd(3).tri=this_mesh.faces;
   % now do check that all shells are contained in each other - dipoli is
   % very picky in this regard...
   % start with brain in skull
   fv=[];
   fv.vertices=bnd(2).pos;
   fv.faces=bnd(2).tri;  
   %in_vert=inpolyhedron(fv,bnd(1).pos);
   % 2020_11 switched to in_polyhedron, inpolyhedron had some random fails
   in_vert=in_polyhedron(fv,bnd(1).pos);
   if sum(~in_vert)>0 
    error('The brain class is not fully contained in the skull (should actually never happen by the way this is generated)')
   end
   
   % skull in scalp
   fv=[];
   fv.vertices=bnd(3).pos;
   fv.faces=bnd(3).tri;  
   in_vert=in_polyhedron(fv,bnd(2).pos);
   if sum(~in_vert)>0 
    warning('The skull class is not fully contained in the scalp - will try to expand the scalp class and try again...')
    if ~isfield(subject,'dipoli_fail')
     subject.dipoli_fail=1;
    elseif subject.dipoli_fail<3
     subject.dipoli_fail=subject.dipoli_fail+1;
    else
     % this is a fail, so delete the modified file to prevent excessive
     % expansion
     fprintf('There is a problem with the headmodel classes not coming to a legal configuration (being fully contained in each other).\n')
     fprintf('Likely there is a problem with the underlying segmentation\n')
     fprintf('Check this with nmri_check_seg (tissue class in particular)\n\n')
     delete(subject.mri_seg);
     % and copy back the original
     movefile(strrep(subject.mri_seg,'.mat','_orig.mat'),subject.mri_seg);
     error('Even after 3 levels of scalp/tissue expansion, the scalp class was still to small for this skull...will stop here and delete the previous attempts to fix by expansion')
    end
    % now expand by 1 voxel
    load(subject.mri_seg,'mri_seg')
    if ~isfield(mri_seg,'dipoli_fail')
     mri_seg.dipoli_fail=1;
    else
     mri_seg.dipoli_fail=subject.dipoli_fail+1;
    end
    scalpdil = imdilate(mri_seg.tissue>0, strel_bol(1)); % expand this a bit
    mri_seg.tissue = scalpdil > 0.5;
    if mri_seg.dipoli_fail==1
     % remember the original file
     movefile(subject.mri_seg,strrep(subject.mri_seg,'.mat','_orig.mat'));
    end
    save(subject.mri_seg,'mri_seg');
    if exist(subject.electrodes_aligned,'file')
     delete(subject.electrodes_aligned)
    end
    % now call nested
    subject=nmri_make_hdm_suma(subject);
   end
   thdm     = bnd;
  case 'simbio'
   % FEM head model
   cfg.method    = 'simbio'; 
   thdm     = fem; 
   cfg.conductivity=zeros(size(fem.tissuelabel));
   % now check the conductivity - Fieldtrip defaults (as of 12/2016)
   for i=1:length(fem.tissuelabel)
    switch fem.tissuelabel{i}
     case 'grey'
      cfg.conductivity(i)=0.33;
     case 'gray' % just in case...
      cfg.conductivity(i)=0.33;
     case 'white'
      cfg.conductivity(i)=0.14; 
     case 'csf'
      cfg.conductivity(i)=1.79; 
     case 'skull'
      cfg.conductivity(i)=0.01; 
     case 'scalp'
      cfg.conductivity(i)=0.43; 
     case 'air'
      cfg.conductivity(i)=0;      
     otherwise
      error(['No conductivity is predefined for class=' fem.tissuelabel{i} ])
    end
   end
  otherwise
   error('unsupported head model for EEG')
  end 
 else
  error('Unsupported datatype (needs to be MEG or EEG currently)') 
 end
 hdm           = ft_prepare_headmodel(cfg,thdm);
end

clear thdm

%% now make leadfield - this should be done always, since we may have different sensors

% get rid of bad channels
cfg           = [];
if isfield(data,'bad_channels')
 good_channels={};
 for i=1:length(data.label)
  if ~any(strcmp(data.label{i},data.bad_channels))
   good_channels(end+1)=data.label(i);
  end
 end
 
 % also filter with elec_present and cross-check with headmodel, mismatch
 % can occur with more recent versions of EGI software having a E257 chanel
 % but this may alos occur otherwise
%  tmp=intersect(good_channels,elec_present.label,'stable');
%  if length(tmp)~=length(good_channels)
%   warning(['There is a mismatch between electrodes present in the headmodel and data. Taking only the common channels.'])
%   good_channels=tmp;
%  end
 
 cfg.channel=good_channels;
 
 % we do not need bad_channels any more now, remove to avoid Fieldtrip
 % warnings
 if isfield(data,'bad_channels')
  data=rmfield(data,'bad_channels');
 end
 
end
data=ft_selectdata(cfg,data);

% leadfield
cfg = [];
if (strcmp(subject.dtype,'MEG') && isfield(data,'grad'))
 cfg.grad    = data.grad;                  % sensor information
elseif (strcmp(subject.dtype,'EEG') && exist('elec_aligned','var'))
 

 cfg.elec    = elec_present;     % electrode positions
else
 error('No gradiometers or electrodes defined')
end
cfg.channel = {'MEG*1'};                %  use {'MEG*2','MEG*3'} for gradiometers
% cfg.channel = 'megmag'; 
% cfg.channel = data.label; 
cfg.sourcemodel.pos = suma_all.pos;              % source points
cfg.headmodel = hdm; 
cfg.normalize = 'no';                    % set to no as of 08/March/2017, recommended by Fieldtrip mailing list
 
leadfield = ft_prepare_leadfield(cfg, data);

%% if we got this far, place a stamp for completed hdm-SUMA
subject=nmri_stamp_subject(subject,['hdmSUMA_' params.headmodel],params);
 
%% save everything incl. subject info
save(subject.hdm_lead,'hdm','leadfield','bnd','subject');

%% visually check the alignement
disp('Generating brain / sensors figure...')
hFig=figure('Position',[0,0,800,800],'Visible','off');
hold on     % plot all objects in one figure
switch params.headmodel
 case 'singleshell'
  % usual MEG model
  ft_plot_headmodel(hdm, 'facecolor','cortex','edgecolor', 'none','facealpha',0.6);
  ft_plot_mesh(bnd(3),'facecolor','skin','edgecolor','none','facealpha',0.2);  
 case 'dipoli'
  ft_plot_mesh(bnd(1),'facecolor','cortex','edgecolor','none','facealpha',0.6);  
  ft_plot_mesh(bnd(2),'facecolor',[1 1 1],'edgecolor','none','facealpha',0.2);  
  ft_plot_mesh(bnd(3),'facecolor','skin','edgecolor','none','facealpha',0.2);  
 case 'simbio'
  ft_plot_mesh(fem, 'surfaceonly', 'yes','facecolor','skin','edgecolor', 'none','facealpha',0.2);
end
% plot the leadfield points
ft_plot_mesh(suma_all.pos,'vertexcolor','b');
% plot fiducials 
nas=[mri_neuromag.cfg.fiducial.nas 1]*mri_neuromag.transform';
lpa=[mri_neuromag.cfg.fiducial.lpa 1]*mri_neuromag.transform';
rpa=[mri_neuromag.cfg.fiducial.rpa 1]*mri_neuromag.transform';
ft_plot_mesh(nas(1:3),'vertexcolor','w','vertexsize',30,'vertexmarker','x');
ft_plot_mesh(lpa(1:3),'vertexcolor',[0 1 0],'vertexsize',30,'vertexmarker','x');
ft_plot_mesh(rpa(1:3),'vertexcolor',[1 0 0],'vertexsize',30,'vertexmarker','x');

if (strcmp(subject.dtype,'MEG') && isfield(data,'grad'))
 ft_plot_sens(data.grad,'style', 'g*');
 clipping=[500 700 150 0];
elseif (strcmp(subject.dtype,'EEG') && exist('elec_aligned','var'))

 ft_plot_mesh(elec_present.chanpos,'vertexmarker','k.','vertexcolor',[0 0 0],'vertexsize',15);
 ft_plot_mesh(elec_missing.chanpos,'vertexmarker','k.','vertexcolor',[1 0 0],'vertexsize',25);
 clipping=[800 800 0 0];
end
hold off
camlight('headlight');

 
% Export figure if wanted
if (isfield(params,'QC_plots') && params.QC_plots==1)
 % make the QC dir
 QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
 if (~exist(QCdir,'dir'))
  mkdir(QCdir)
 end
 
 hFig=nmri_rotate_print(hFig,[0 0;180 0;90 0;0 90],fullfile(QCdir,['hdm_sensors_plot_' subject.id '_' subject.exam_id '_' params.headmodel]),[subject.id ' - Overview'],clipping);
end
 
% only show after save
% title(subject.id);
% view([0 0])
% set(hFig,'Visible','on')
 
disp('...Headmodelling done')
disp('Check the alignement now and proceed to nmri_processing when happy'); 
 
end

