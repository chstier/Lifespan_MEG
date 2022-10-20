function [ gta_metrics ] = nmri_calculate_gta( full_metrics, opt)
% [ gta_metrics ] = nmri_calculate_gta( full_metrics, opt)
%   Will calculate GTA metrics (using BCT functions) 
%   allows normalization by permutating
%   
%
% full_metrics      = cell array of metrics, e.g. from nmri_export_metrics
%
%
% opt               = struct array of options
%  .metric          = GTA metric to calculate (using BCT functions)
%                     'cc_wu': clustering_coef_wu
%                     'cc_bu': clustering_coef_bu
%                     'cc_bd': clustering_coef_bd
%                     'cc_wd': clustering_coef_wd
%                     'ec': eigenvector_centrality_und
%                     'cpl_wei': distance_wei + charpath
%
% .normalize        = true/false: do normalization by permutation
% .normalize_perm   = number of random permutations, default: 100
%
% gta_metrics  = the calcuated GTA metrics (cell array)
 

% check opt and set defaults
if ~exist('opt','var')
 % defaults
 opt=[];
end
if ~isfield(opt,'metric') || ~ischar(opt.metric) || isempty(opt.metric)
 error('need to provide a metric to calculate')
end

if ~exist(opt.metric,'file')==2
 error(['Function = ' opt.metric ' not found. Check the path to BCT (needs to be in Matlab path).'])
end

if ~isfield(opt,'normalize_perm') || ~isnumeric(opt.normalize_perm) || opt.normalize_perm<1 
 opt.normalize_perm=100;
end

if ~isfield(opt,'normalize') || ~islogical(opt.normalize)  
 opt.normalize=true;
end

N=length(full_metrics);
gta_metrics=cell(N,1);
fprintf('Metric: %s\n',opt.metric)

for n=1:N
 % now loop per subjects
 this_full= full_metrics{n};
 % remove NaN
 this_full(isnan(this_full))=0;
 fprintf('Subject: %d/%d\n',n,N)
 
 % calculate
 gta_metrics{n}=calc_gta(this_full,opt.metric);
 
 % normalize?
 if opt.normalize
  ref_norm=zeros(size(gta_metrics{n}));
  fprintf('Normalizing (%d permutations)...\n',opt.normalize_perm)
  for i=1:opt.normalize_perm
   shuffeled=calc_gta(sub_shuffleAdjacencyMatrix(this_full),opt.metric);
   fprintf('.')
   ref_norm=ref_norm+shuffeled;
  end
  ref_norm=ref_norm/opt.normalize_perm;
  gta_metrics{n}=gta_metrics{n}./ref_norm;
  fprintf('done\nRescaled GTA metric to mean=%f\n',mean(gta_metrics{n}))
 end

end
end

% local functions
function gta_metric=calc_gta(full_metric,metric)
 switch metric
  case 'cc_wu'
   gta_metric=clustering_coef_wu(full_metric);
  case 'cc_wd'
   gta_metric=clustering_coef_wd(full_metric);
  case 'cc_bu'
   gta_metric=clustering_coef_bu(full_metric);
  case 'cc_bd'
   gta_metric=clustering_coef_bd(full_metric);
  case 'ec'
   gta_metric=eigenvector_centrality_und(full_metric);
  case 'cpl_wei'
   gta_metric=charpath(distance_wei(full_metric));

  otherwise
    error(['Unknown metric: ' metric])
 end
end
