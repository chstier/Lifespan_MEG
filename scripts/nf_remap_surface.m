function [ out_vertices, out_faces, vertex_mapping ] = nf_remap_surface( in_vertices, in_faces )
%[ out_vertices, out_faces, vertex_mapping ] = nf_remap_surface( in_vertices, in_faces )
% remove vertices set to NaN (at y coord) from surface and relabels face
% indices

% method from http://www.gomatlab.de/vertices-eckpunkte-mit-nan-aus-faces-flaechen-entfernen-t36528.html

% Select vertices to be deleted (marked with NaN) 
remove = any(isnan(in_vertices(:,2)),3); 
keep = ~remove; 
out_vertices = in_vertices(keep, :);


% Remove faces, which have any deleted point: 
removeInd = find(remove); 
newFaces = in_faces(~any(ismember(in_faces, removeInd), 2), :); 

% Clean up the indices of faces: 
keepInd = cumsum(keep); 
out_faces = keepInd(newFaces); 

% generate vertex to vertex mapping
vertex_mapping=[1:length(in_vertices)]';
vertex_mapping(remove)=[];

end

