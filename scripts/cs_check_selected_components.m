function [subject, data] = cs_check_selected_components(subject, data)
% This function loads the selected ICA components for each channel type for
% a visual check
% written by Christina

% check the call
if (~exist('subject','var') ) 
 error('Need a valid subject struct or .m/.mat file to work with')
end

% call the subject and params include
nmri_include_read_ps

%% get data of selected components for each channel type 
chan = {'mag', 'grad'};

for c = 1:length(chan)
 datname = fullfile([chan{c} '_selected_comp_data']);
 subject.(datname)=fullfile(subject.analysis_dir,subject.id,'processed',['ICA_sel_comp_' chan{c} '_' subject.id '_' subject.exam_id '.mat']);
 comp_data.(chan{c}) = load(subject.(datname), 'data'); % store in struct
end

%% make the QC dir
if (~isfield(subject,'QCdir'))
 subject.QCdir=fullfile(subject.analysis_dir,subject.id,'QC');
end
if (~exist(subject.QCdir,'dir'))
 mkdir(subject.QCdir)
end

%% do plotting for review
% MAG components
cfg          = [];
cfg.viewmode = 'component'; 
cfg.layout   = 'neuromag306mag_helmet.mat'; % specify the layout file that should be used for plotting
f = ft_databrowser(cfg, comp_data.mag.data);


% GRAD components
cfg          = [];
cfg.viewmode = 'component';
cfg.layout   = 'neuromag306planar_helmet.mat'; % specify the layout file that should be used for plotting
f2 = ft_databrowser(cfg, comp_data.grad.data);



