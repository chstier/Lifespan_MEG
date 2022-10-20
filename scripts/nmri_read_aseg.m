function [ surf ] = nf_read_aseg( fspath, conf_file )
%[ surf ] = nf_read_aseg( fsid, conf_file )
%   Will read in the aseg / subcortical nuclei as meshes
%   
%   fspath    = path to Freesurfer dir
%   conf_file = path to conf_file (default,
%   '/tools/common/freesurfer-labels-aseg-connectome-10')
%
% NOTE: the meshes are in VOXEL space, i.e. as read in the image
%
%   will return a struct with
% .pos = vertices
% .tri = faces
% .annot = annotations
% .trans_aseg = transformation read by load_mgh from aseg image. Note, this
% may be different from T1.nii. You have to deal with that seperately
%
% NF 11/2016

if ~exist('fspath','var')
 error('Need Freesurfer dir')
end

if ~exist('conf_file','var')
 conf_file='/tools/common/freesurfer-labels-aseg-connectome-10';
end

if ~exist(conf_file,'file')
 error(['Could not find conf file = ' conf_file])
end

aseg_file=fullfile(fspath,'mri','aseg.mgz');
if ~exist(aseg_file,'file')
 error(['Could not find aseg image = ' aseg_file])
end

% load aseg
[aseg_vol aseg_trans] = load_mgh(aseg_file);
% strangly enough the transformation matrix returend by load_mgh does not
% seem to be correct. There is a linear shift compared to the T mat loaded
% from T1.nii (SUMA). Not clear, which one is really correct. So, for now
% no transformation is applied here. It will be up to the super-function to
% deal with this. The meshes are, thus, returned in voxel-space 

% load our configuration file - IDval, name, Nvert, hemi, erode in vx
fid=fopen(conf_file,'r');
aseg=textscan(fid,'%d %s %d %d %d %*[^\n]','CommentStyle','#','MultipleDelimsAsOne',1,'Delimiter',' ');
fclose(fid);

% check the sizes
if (length(aseg)<5)
 conf_file
 error('Configuation file did not give the neceassary info, check')
end


surf=[];
surf.pos=[];
surf.tri=[];
surf.annot=[];
surf.annot_key=cell(1,2);
surf.annot_key{1}=zeros(size(aseg{1}));
surf.annot_key{2}=cell(size(aseg{1}));
surf.cortex=[];
surf.hemi=[];
surf.msk=[];
surf.trans_aseg=vox2ras_0to1(aseg_trans); % safe as 1-based transform
% now we shall parse all ROIs
for i=1:length(aseg{1})
 % get the current ROI
 this_roi=double(aseg_vol==aseg{1}(i)); 
 % erode image somewhat
 erK=nf_strel_bol(double(aseg{5}(i)));
 this_roi=imerode(this_roi,erK);
 this_roi=nf_volumesmooth(this_roi,2);
 % rethreshold
 %this_roi=double(this_roi>0.6);
 % make a mesh
 [this_mesh.faces, this_mesh.vertices] = isosurface(this_roi, 0.6);
 this_mesh.vertices = this_mesh.vertices(:,[2 1 3]); % matlab is twisted...
 
 % option for iso2mesh, seems to work much slower and less nice
% opt = [];
% opt.radbound = 3; % set the target surface mesh element bounding sphere be <3 pixels in radius
% opt.maxnode = aseg{3}(i);
% opt.dofix = 1;
% [pos, tri] = v2s(this_roi, 1, opt, 'cgalsurf');
% tri = tri(:,1:3);
% this_mesh.vertices = pos;
% this_mesh.faces = tri;
 
 % smooth with iso2mesh function
 this_mesh.vertices=smoothsurf(this_mesh.vertices,[],meshconn(this_mesh.faces,length(this_mesh.vertices)),10,1,'lowpass');
 
 % decimate to target
 target_vertices=aseg{3}(i);
 target_faces=target_vertices*2; % some first target
 while(length(this_mesh.vertices)>target_vertices)
  [this_mesh.faces, this_mesh.vertices]=reducepatch(this_mesh,target_faces);
  % reduce by further 2% if needed
  target_faces=floor(length(this_mesh.faces)*0.98);
  fprintf('Reducing mesh (vx count=%d)\n',length(this_mesh.vertices))
 end
 % warp the mesh
 % not implented here

 % print some info
 fprintf('For ROI=%s we have %d vertices\n\n',aseg{2}{i},length(this_mesh.vertices));
 
 % now add to summary mesh and make .pos / .tri (Fieldtrip-style)
 oldVc=length(surf.pos);
 surf.pos=vertcat(surf.pos,this_mesh.vertices);
 surf.tri=vertcat(surf.tri,(this_mesh.faces+oldVc));
 surf.annot=vertcat(surf.annot,repmat(aseg{1}(i),length(this_mesh.vertices),1)); % keep Freesurfer key as annot
 surf.annot_key{1}(i)=aseg{1}(i);
 surf.annot_key{2}{i}=aseg{2}{i};
 surf.cortex=vertcat(surf.cortex,repmat(2,length(this_mesh.vertices),1)); % all non-cortex
 surf.msk=vertcat(surf.msk,repmat(1,length(this_mesh.vertices),1)); % all in mask
 surf.hemi=vertcat(surf.hemi,repmat(aseg{4}(i),length(this_mesh.vertices),1)); % keep hemi
end



end

