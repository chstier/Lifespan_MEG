function nf_json_write( filename, data )
%nf_json_write( filename, data )
%   Writes a JSON data file with some formatting to be editable
 fid=fopen(filename,'w');
 [json, delim]=split(jsonencode(data),{',','{','(','}',']','":'});
 charcount=0;
 for i=1:length(json)-1 % last item should always be empty
  if (strcmp(delim{i},'":'))
   fprintf(fid,'\n');
   charcount=0;
  end
  fprintf(fid,'%s',json{i});
  if strcmp(delim{i},'}')
   fprintf(fid,'\n');
   charcount=0;
  end
  fprintf(fid,'%s',delim{i});
  charcount=charcount+length(json{i});
  %if i<length(json) && ( charcount>80 || strcmp(delim{i},'}') || strcmp(delim{i},']'))
  % fprintf(fid,'\n');
  % charcount=0;
 % end
 end
 fclose(fid);
end

