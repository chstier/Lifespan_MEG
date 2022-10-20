function basename = get_basename(fullfilename)
% returns basename from filename (if standard convention P/C.name.seq)
% written by Niels Focke 04/2009 for new UMG basename code

basename='';

% check if we have / in fullfilename
slashes=strfind(fullfilename,'/');

if (~isempty(slashes))
 [pa, fi, ext] = fileparts(fullfilename);
 if isempty(fi) 
  error('filename appears bad...could not extract file from full path');
 end
else
 fi=fullfilename;
end



k = strfind(fi, 'D.ABCD01');
if (~isempty(k))
 basename='D.ABCD01';
else
 % try BIDS style first
 if ~isempty(regexp(fi,'(sub-.*?)(_ses-.*?){0,1}_'))
  basename=char(regexp(fi,'(sub-.*?)(_ses-.*?){0,1}_','match'));
  % remove trailing _
  basename=basename(1:end-1);
 else
  % use storage logic
  basename=char(regexp(fi,'[CPSO]\.[A-Z][A-Z][A-Z][A-Z][0-9][0-9][A-Z]{0,1}','match'));
 end
end

