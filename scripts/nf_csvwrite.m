function nf_csvwrite( file, data )
%nf_csvwrite( file, data )
%   Dumps 2D cell array (mixed possible) into a CSV file

if (~iscell(data) || ndims(data)~=2)
 error('Data needs to be a 2D cell array')
end

af=fopen(file,'w');

for row=1:size(data,1)
 for col=1:size(data,2)
  if (ischar(data{row,col}))
   fprintf(af,'%s',data{row,col});
  else
   % seems to be number
   number=data{row,col};
   if (number==floor(number))
    % integer number
    fprintf(af,'%g',number);
   elseif abs(number)<1e-6
    % print with exp for very small
    fprintf(af,'%e',number);
   else
    % print with .
    fprintf(af,'%f',number);
   end
  end
  if (col<size(data,2))
   fprintf(af,'\t');
  end
 end
 fprintf(af,'\n');
end

fclose(af);
end

