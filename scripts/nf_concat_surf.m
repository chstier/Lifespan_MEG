function [ surf_con ] = nf_concat_surf( surf_A, surf_B )
%[ surf_con ] = nf_concat_surf( surf_A, surf_B )
%   Vertically concatenates two surface struct (taking care of special fields)
%   expects .faces/.vertices (Matlab-style) or .pos/.tri (Fieldtrip style)

f_A=fieldnames(surf_A);
f_B=fieldnames(surf_B);

f_all=union(f_A,f_B);

%if (length(f_A)~=length(f_B))
% error('Your surfaces to not have the same number of struct fields...this will not work')
%end

surf_con=[];
for i=1:length(f_all)
 this_f=f_all{i};
 if isfield(surf_A,this_f) 
  this_A=surf_A.(this_f);
 else
  warning(['Could not find field=' this_f ' in surface_A - skipping this field'])
  continue
 end
 if isfield(surf_B,this_f) 
  this_B=surf_B.(this_f);
 else
  warning(['Could not find field=' this_f ' in surface_B - skipping this field'])
  continue
 end
 
 if strcmpi(this_f,'tri') || strcmpi(this_f,'vertices') 
  % faces are special, we need to shift by length of A
  if (isfield(surf_A,'vertices'))
   shift_by=length(surf_A.vertices);
  elseif (isfield(surf_A,'pos'))
   shift_by=length(surf_A.pos);
  else
   error('Could not find a .vertices or .pos field -- cannot shift faces/triangles')
  end
  % now shift faces
  this_B=this_B+shift_by;
 end

 % now concatenate
 if strcmpi(this_f,'annot_key')
  % special for annot key  
  surf_con.annot_key=cell(1,2);
  for ii=1:2
   surf_con.annot_key{ii}=vertcat(surf_A.annot_key{ii},surf_B.annot_key{ii});
  end
 else
  % default
  surf_con.(this_f)=vertcat(this_A,this_B);
 end
end

end


