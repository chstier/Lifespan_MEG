function [ mff_reader ] = nmri_check_mff_reader( mff )
% [ mff_reader ] = nmri_check_mff_reader( mff )
%   Will auto-determine the MFF reader to use

if ~exist(mff,'file')
 error('Could not find mff=%s',mff)
end

xmlfile=fullfile(mff,'info.xml');
if ~exist(xmlfile,'file')
 error('Could not find info.xml=%s',xmlfile)
end

xml = nf_xml2struct(xmlfile);

if isfield(xml,'mffVersion') && str2double(xml.mffVersion)>2
 % at least version 3 file, so use JAR reader
 mff_reader='egi_mff_v3'; % changed to V3
else
 % legacy file (<ver 3), need to use the old v1 reader
 mff_reader='egi_mff_v1';
end

host=getenv('HOST');
if strcmp(host,'nmri-srv.be-mrz.klin')
 % for memory reasons, always use _v1 on NMRI server, too high mem  load
 % with _v3 reader
 mff_reader='egi_mff_v1';
end

end

