function [ vertices, faces ] = nf_load_surf_flexibele( surf_file )
%[ vertices, faces ] = nf_load_surf_flexibale( surf_file )
% Wrapper script to load the different surfaces formats, i.e.
% Freesurfer format (no extension)
% ASCII (.asc extension)
% GIFTI (.gii extension)

%check if we have the files
if (~exist(surf_file,'file'))
 error(['Surface file (= ' surf_file ') not found'])
end

[~, ~, ext]=fileparts(surf_file);

ext=lower(ext); % make lower case

if (strcmp(ext,'.asc') || strcmp(ext,'.ascii'))
 % ASCII format assumed
 disp('Using ASCII reader')
 [vertices, faces]=freesurfer_read_ascii(surf_file);
elseif (strcmp(ext,'.gii') || strcmp(ext,'.gifti'))
 % GIFTI, use reader included in SPM
 disp('Using SPM-GIFTI reader')
 g=gifti(surf_file);
 vertices=g.vertices;
 faces=g.faces;
else
 % assume default Freesurfer format
 disp('Using Freesurfer reader')
 [vertices, faces]=freesurfer_read_surf(surf_file);
end


end

