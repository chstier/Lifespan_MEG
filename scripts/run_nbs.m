function run_nbs(filename,threshold,designMatrix,nodecoor,nodelabels)


h = figure;
close(h)
clear nbs
clear UI

UI.alpha.ui    = '0.05';
% Scalar
% Significance (alpha threshold).
% See also NBSstats

UI.method.ui   = 'Run NBS'; %| 'Run FDR'
% Perform NBS or FDR?


UI.test.ui     = 't-test'; %| 'F-test' |  'One Sample'
% Statistical test to perform.
% See also NBSglm

UI.size.ui     = 'Extent'; %'Extent'; %| 'Intensity'
% Use intensity or extent to measure size of a
% network component?
% [optional if UI.method.ui=='Run FDR']
% See also NBSstats


UI.perms.ui    = '5000';
% Scalar integer
% Number of permutations.
% See also NBSglm



UI.contrast.ui = '[-1 1]';
% 1 x p numeric array specifying contrast, where p
% is the number of independent variables in the GLM.
% Must be specified as a valid Matlab expression
% for a 1 x p array
% See also NBSglm


UI.design.ui   = designMatrix;

% n x p numeric array specifying a design matrix,
% including a column of ones if necessary. p is the
% number of independent variables in the GLM, n is
% the number of observations.
% Can be specified either as a:
%1. Valid Matlab expression for an n x p array
%2. Text file containing numeric data arranged into
%   n rows and p columns
%3. A binary Matlab file (.mat) storing an n x p
%   numeric array
%See also NBSglm


UI.matrices.ui =   filename;   
% N x N numeric array specifying a symmetric
% connectivity matrix for each of M observations
%(e.g. subjects), where N is the number of nodes.
% Can be specified either as a:
%1. Valid Matlab expression for an N x N x M array
%2. A total of M seperate text files stored in a
%   common directory, where each text file contains
%   numeric data arranged into N rows and N columns.
%   Specify only one such text file and the others
%  within the same directory will be identified
%  automatically.
%3. A binary Matlab file (.mat) storing an N x N x M
%   numeric array


UI.exchange.ui   = '';
UI.node_coor.ui  = nodecoor;
UI.node_label.ui = nodelabels;


UI.thresh.ui =  threshold;
% Primary test statistic threshold.
%[optional if UI.method.ui=='Run FDR']
% See also NBSstats

NBSrun(UI,h)



end