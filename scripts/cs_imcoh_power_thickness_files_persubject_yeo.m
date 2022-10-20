%% This script loads the PALM group files and vertex-wise connectivity, power, and cortical thickness values, does resampling for the yeo-resolution and saves
% .mat-files in format 
% 
% Written and modified by Christina Stier, 2022

%% create files per frequency band and and metric - regional values 
% nSubj x regimcoh(16x1) x regpow(16x1) x reg

%load SUMA
load('conf/atlas/Yeo2011_7Networks_N1000_suma-all-fsaverage-10.mat','suma_all');

metric = {'coh_img_', 'power_', 'thickness_smoothed_12_'};
freq = {'Delta_','Theta_','Alpha_','Beta1_','Beta2_','Gamma_'};
scale = {'abs_not_scaled_N2'};
rootdir_MEG = ('/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohortX/'); % take all subjects
rootdir_th = ('/home/uni10/nmri/projects/cstier/aging_pipeline/export/cohortX/October-2022_ld10_sm12/');
outputname = 'subjmatrix_I_P_T_yeo_';
group = 'all';


for var_freq = 1:length(freq)
   % load values for imcoh
   for var_metric =  metric{1}
       data = load([rootdir_MEG...
            metric{1}...
            freq{1, var_freq}...
            scale{1} '.mat']);
        for i = 1:length(data.scale_metrics)
            m1 = data.scale_metrics{i};
            m1(2005:2338,1) = NaN;
            regvalues_imcoh = cs_nf_ft_report_regional(suma_all,m1);
            regvalues_imcoh(1,:) = []; 
            roi = regvalues_imcoh(:,1); % get the ROIs to append later
            regvalues_imcoh(:,1:3) = [];
            subj(i).regimcoh = regvalues_imcoh;
            %regval{:,1,i} = regvalues_imcoh;
        end
   end
   
    % load values for power
    for var_metric =  metric{2}
            data = load([rootdir_MEG...
            metric{2}...
            freq{1, var_freq}...
            scale{1} '.mat']);
        for i = 1:length(data.scale_metrics)
            m2 = data.scale_metrics{i};
            m2(2005:2338,1) = NaN;
            regvalues_power = cs_nf_ft_report_regional(suma_all,m2);
            regvalues_power(1,:) = []; 
            roi = regvalues_power(:,1);
            regvalues_power(:,1:3) = [];
            subj(i).regpow = regvalues_power;
            %val{1,2,i} = regvalues_power;
        end
    end
    
    % load thickness
    for var_metric =  metric{3}
            data = load([rootdir_th...
            metric{3}...
            'suma_all_N2.mat']);
        for i = 1:length(data.thickness) % change here...
            m3 = data.thickness{i};
            m3(2005:2338,1) = NaN;
            regvalues_th = cs_nf_ft_report_regional(suma_all,m3);
            regvalues_th(1,:) = []; 
            roi = regvalues_th(:,1);
            regvalues_th(:,1:3) = [];
            subj(i).regth = regvalues_th;
        end
    end
    
   filename = strcat(outputname, cell2mat(freq(1,var_freq)), group);
   save(filename, 'subj')
  
end

