function nf_coreg_single(modus,targetfile,referencefile,outdir)
% function nf_coreg_single(modus,targetfile,referencefile,outdir)
% coregs one image to one reference image
% modus         = 1: only estimation, information will be stored in header and new image will be in outdir
%                 2: estimate and reslice, both files will be in outdir (trilinear)
%                 3: estimate and reslice (nearest)
%                 4: estimate and reslice (B-spline)
%                 5: resclice only (trilinear)
%                 6: resclice only (nearest)
%                 7: resclice only (B-spline)
% targetfile 	= a single file
% referencefile	= a single reference image
% outdir    	= output dir (1st step is to copy source files here)
%
% written by Niels Focke 11/2016


global defaults last_root_dir

if ( ~isfield(defaults,'modality'))
 spm('PET');
end


fh=spm_figure('FindWin','Interactive');

flags=defaults.coreg;

if ( ~exist('modus') || ~isnumeric(modus) || modus<1 || modus>5 )
 modus=spm_input('Coreg mode','1','m','Estimate (and store in header)|Estimate and reslice (trilinear)|Estimate and reslice (nearest)|Estimate and reslice (B-spline)',[1 2 3 4]);
end

if ( ~exist('targetfile') || isempty(targetfile) )
 targetfile=spm_select(1, 'image', 'Choose target image (the one to be coregistered)');
end

tvols = spm_vol(targetfile);
N = length(tvols);

if ( ~exist('referencefile') || isempty(referencefile) )
 referencefile=spm_select(1, 'image', 'Choose reference image (the one that will stay the same)');
end

if ( ~exist('outdir') || isempty(outdir) )
 outdir=spm_select(1, 'dir', 'Choose destination dir (where output should be written to)');
end
if ( isempty(outdir) )
 outdir=[cd '/'];
end

for n = 1:N
 bn=get_basename_or_root(char(tvols(n).fname));

 
 if (modus==1)
  [pa, fi, ext] = fileparts(tvols(n).fname);
  savefile=strcat(pa,'/',fi,'.*');
  system(['cp -L ' savefile ' ' outdir]);
  tvols(n).fname=[outdir fi ext];
 end
 
 switch modus
 case 1
  matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.estimate.eoptions=flags.estimate;
 case 2
  matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions=flags.estimate;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions=flags.write;
 case 3
  matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions=flags.estimate;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions=flags.write;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp=0;
 case 4
  matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions=flags.estimate;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions=flags.write;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp=2;
 case 5
  matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.write.roptions=flags.write;
 case 6
  matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.write.roptions=flags.write;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp=0;
 case 7
  matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(referencefile);
  matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(tvols(n).fname);
  matlabbatch{1}.spm.spatial.coreg.write.roptions=flags.write;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp=2;
 end
 spm_jobman('run',matlabbatch);

 if (modus~=1)
  [pa, fi, ext] = fileparts(tvols(n).fname);
  savefile=strcat(pa,'/r',fi,'.*');
  if ( ~strcmp(pa,outdir) ) && ( ~strcmp([pa '/'],outdir) ) 
   movefile(savefile,outdir); 
  end
 end
 
end


