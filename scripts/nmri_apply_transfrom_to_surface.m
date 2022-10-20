function [ vertices_out ] = nmri_apply_transfrom_to_surface( vertices_in , y_file )
%[ surf_out ] = nmri_apply_transfrom_to_surface( surf, y_file )
%   Applies a SPM12 y_ transform (point-point mapping) to a surface in
%   Fieldtrip or Matlab style (will transform verticies / points)


if ~exist('y_file','var') || ~exist(y_file,'file')
 error('No y_file specified')
end

% require SPM12 to be present
defaults=nf_request_spm;

vv=[[y_file ',1,1']; [y_file ',1,2']; [y_file ',1,3']];
vols=spm_vol(vv);
[mni, XYZ]=spm_read_vols(vols);

vertices_out=zeros(size(vertices_in));
for i=1:size(vertices_in,1)
 x=round(vertices_in(i,1))+1;
 y=round(vertices_in(i,2))+1;
 z=round(vertices_in(i,3))+1;

 vertices_out(i,1)=mni(x,y,z,1);
 vertices_out(i,2)=mni(x,y,z,2);
 vertices_out(i,3)=mni(x,y,z,3);
end

end

