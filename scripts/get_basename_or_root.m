function bn = get_basename_or_root(fullfilename)
% returns basename from filename (if standard convention P/C.name.seq)
% if not possible try to determine root filename an return this
% written by Niels Focke 05/2009 for new UMG basename code

bn=get_basename(char(fullfilename));
if (strcmp(bn,'')) 
 % if no basename can be found use full filename
 [pa, bn, ext] = fileparts(strtok(char(fullfilename),','));
  %shaving off known SPM prefixes to estimate orginal filename
  shave=true;
  while(shave)
   if (strncmp(bn,'c1',2))
    bn=bn(3:end);
   elseif (strncmp(bn,'c2',2))
    bn=bn(3:end);
   elseif (strncmp(bn,'c3',2))
    bn=bn(3:end);
   elseif (strncmp(bn,'x',1))
    bn=bn(2:end);
   elseif (strncmp(bn,'w',1))
    bn=bn(2:end);
   elseif (strncmp(bn,'m',1))
    bn=bn(2:end);
   elseif (strncmp(bn,'s',1))
    bn=bn(2:end);
   elseif (strncmp(bn,'r',1))
    bn=bn(2:end);
   elseif (strncmp(bn,'i',1))
    bn=bn(2:end);
   elseif (strncmp(bn,'JM',2))
    bn=bn(3:end);
   else
    shave=false;
   end 
  end
  
  % it may happen that we have _seg_sn at the end
  pos=strfind(bn,'_seg_sn');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  
  pos=strfind(bn,'_seg_inv_sn');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  
  pos=strfind(bn,'_sn');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  
  % shave off image tags, lasy fix but should work in most cases
  pos=strfind(bn,'.cFLAIR');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  pos=strfind(bn,'.FLAIR');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  pos=strfind(bn,'.hrT1');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  pos=strfind(bn,'.DTI');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  pos=strfind(bn,'.EPI_DTI');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  pos=strfind(bn,'.T2map');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  pos=strfind(bn,'.FDG_PET');
  if (length(pos)>0)
   bn=bn(1:pos-1);
  end
  
end

