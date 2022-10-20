function [ hFig ] = nmri_schemaball_plot( metric, ref_surf, opt )
%[ hFig ] = nmri_schemaball_plot( metric, ref_surf, opt )
%   will generate a schemaball plot from a metric, with some options
% opt.scale = number to scale the metric from -1 to 1, i.e. the maximal
%             value. If not provided, will take abs-max of metric
% opt.show_labels = 0: hide labels, 1: show full labels, 2: abbr. labels
% opt.node_scale = 0: no scaling of nodes, 1: do scale (default)
% opt.ccolor: link color (optional)
% opt.ncolor: node color (optional)


% check opt and set defaults
if ~exist('opt','var')
 % defaults
 opt=[];
end
if ~isfield(opt,'scale') || ~isnumeric(opt.scale) || isempty(opt.scale)
 opt.scale=max(abs(metric(:)));
end
if ~isfield(opt,'show_labels') || ~isnumeric(opt.show_labels) || isempty(opt.show_labels)
 opt.show_labels=1;
end
if ~isfield(opt,'node_scale') || ~isnumeric(opt.node_scale) || isempty(opt.node_scale)
 opt.node_scale=1;
end

%if ~isfield(opt,'lcolor') || ~isnumeric(opt.lcolor) || isempty(opt.lcolor)
% opt.node_scale=1;
%end

% check surface
if ~exist('ref_surf','var') || ~isstruct(ref_surf) || ~isfield(ref_surf,'annot') || ~isfield(ref_surf,'annot_key') 
 error('The reference surface does not have the required format')
end


% check metric
if ~exist('metric','var')
 error('You need to provide the full brain metric')
elseif ~isnumeric(metric)
 error('Metric needs to be given as numeric array (single subject)')
end


% assume that we need to flip around the 2nd half of the metric to 
N=length(metric);
if mod(N,2)~=0
 error('Metric is not even, seems to be asymmetric - cannot flip around the 2nd half')
end
N2=N/2;

% flip lower half
metric(N2+1:end,:)=flip(metric(N2+1:end,:),1);
metric(:,N2+1:end,:)=flip(metric(:,N2+1:end),2);

% flip annotion
ref_surf.annot_key{1}(N2+1:end)=flip(ref_surf.annot_key{1}(N2+1:end));
ref_surf.annot_key{2}(N2+1:end)=flip(ref_surf.annot_key{2}(N2+1:end));
ref_surf.annot_key{3}(N2+1:end)=flip(ref_surf.annot_key{3}(N2+1:end));

% deal with scale
metric(metric>opt.scale)=opt.scale;
metric(-metric>opt.scale)=-opt.scale;


% deal with labels
switch opt.show_labels
 case 1
  labels=ref_surf.annot_key{2};
 case 2
  % abbreviate
  labels=ref_surf.annot_key{2};
  wds=regexp(labels,'\w+','match');
  for i=1:length(wds)
   this_l='';
   for ii=1:length(wds{i})
    this_l=[this_l wds{i}{ii}(1:min(2,length(wds{i}{ii})))];
   end
   labels{i}=this_l;
  end
 
 otherwise
 labels=repmat({''},length(ref_surf.annot_key{2}),1);
end

% deal with link color
if ~isfield(opt,'ccolor')
 opt.ccolor=[];
end
                     
% deal with node color
if ~isfield(opt,'ncolor')
 % auto detect 
 if size(ref_surf.annot_key,2)>2
  % we seem to have lobes info
  opt.ncolor=[];
  for i=1:length(ref_surf.annot_key{3})
   switch ref_surf.annot_key{3}{i}
    case 'F'
     this_c=[1 1 0];
    case 'T'
     this_c=[0.5 0.5 1];
    case 'P'
     this_c=[0 1 0.5];
    case 'O'
     this_c=[0.5 1 0];
    case 'I'
     this_c=[0 1 0];
    case 'C'
     this_c=[1 0 1];
    case 'B'
     this_c=[1 0.5 0.5];
    otherwise
     this_c=[1 0 0];
   end
   opt.ncolor=[opt.ncolor; this_c]; 
  end
  
 else
  opt.ncolor = [0 0 1];   %all blue
 end
end

% deal with text color
if ~isfield(opt,'tcolor')
 opt.tcolor=opt.ncolor;
end

hFig=schemaball(metric/opt.scale,labels,opt.ccolor,opt.ncolor,opt.tcolor,opt.node_scale);

end

