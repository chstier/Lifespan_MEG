function [smoothed_right, smoothed_left] = cs_smooth_thickness_suma_all(subject, opt, suma_folder, mod)
%% This function loads individual freesurfer-files (thickness.niml-dset for each hemisphere)
% and smoothes individual maps using AFNI https://afni.nimh.nih.gov/

% Written and modified by Christina Stier, 2020/22

% get smoothing info
sigma_value = opt.sigma;
sm_param = opt.fwhm;
suma_dir = fullfile(suma_folder, [subject.id mod], 'SUMA');
subject.id2 = subject.id;
ld = opt.ld; % choose density factor

%% Smooth thickness-maps

% left hemisphere
if exist(fullfile(suma_dir,['std.' ld '.lh.thickness.niml.dset']), 'file')
 if ~exist(fullfile(suma_dir,['std.' ld '.lh.thickness_smh7_fwhm_' sm_param '.niml.dset']), 'file')
   command = ['SurfSmooth -met HEAT_07 -spec std.' ld '.' subject.id2 mod '_lh.spec -surf_A lh.smoothwm.gii -input std.' ld '.lh.thickness.niml.dset -fwhm ' sm_param ' -sigma ' sigma_value ' -output std.' ld '.lh.thickness_smh7_fwhm_' sm_param]
   [cmdout] = system(['cd ' suma_dir ';' command])
   if exist(fullfile(suma_dir,['std.' ld '.lh.thickness_smh7_fwhm_' sm_param '.niml.dset']), 'file')
    disp('smoothing successful')
   else
    error('smoothing failed')
   end
 else
  disp(['lh smoothing already present for ' subject.id])
 end
else
 error(['lh.thickness.niml.dset file not found for ' subject.id])
end

% right hemisphere
if exist(fullfile(suma_dir,['std.' ld '.rh.thickness.niml.dset']), 'file')
 if ~exist(fullfile(suma_dir,['std.' ld '.rh.thickness_smh7_fwhm_' sm_param '.niml.dset']), 'file')
   command = ['SurfSmooth -met HEAT_07 -spec std.' ld '.' subject.id2 mod '_rh.spec -surf_A rh.smoothwm.gii -input std.' ld '.rh.thickness.niml.dset -fwhm ' sm_param ' -sigma ' sigma_value ' -output std.' ld '.rh.thickness_smh7_fwhm_' sm_param]
   [status cmdout] = system(['cd ' suma_dir ';' command])
   if exist(fullfile(suma_dir,['std.' ld '.rh.thickness_smh7_fwhm_' sm_param '.niml.dset']), 'file')
    disp('smoothing successful')
   else
    error('smoothing failed')
   end
 else
  disp(['rh smoothing already present for ' subject.id])
 end
else
 error(['rh.thickness.niml.dset file not found for ' subject.id])
end

%% name the output-files
smoothed_right = fullfile(['std.' ld '.rh.thickness_smh7_fwhm_' sm_param '.niml.dset']);
smoothed_left = fullfile(['std.' ld '.lh.thickness_smh7_fwhm_' sm_param '.niml.dset']);

end