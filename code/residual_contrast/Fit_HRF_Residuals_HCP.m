function [Results, Flags] = Fit_HRF_Residuals_HCP(Task, Direction)


addpath(genpath('/dcl01/smart/data/CANlab_tools/CanlabCore-master'))      
addpath(genpath('/dcl01/smart/data/CANlab_tools/spm12-master')) 
addpath('/dcl01/smart/data/CANlab_tools/Neuroimaging_Pattern_Masks-master/Atlases_and_parcellations/2018_Wager_combined_atlas')
addpath('/users/mlindqui/HRFresults')
addpath('/users/mlindqui/HRF_Est_Toolbox4')


TR = 0.72;  
T = 30;

%Task = 'WM';
%Task = 'LANGUAGE';
%Direction = 'LR';
%Direction = 'RL';

Type = 'fsl_ar';
%Type = 'fsl_noar';


switch Task

    case 'EMOTION' 
        numstim = 2;
    case 'GAMBLING' 
        numstim = 3;
    case 'LANGUAGE' 
        numstim = 6;
    case 'MOTOR'
        numstim = 6;
    case 'RELATIONAL' 
        numstim = 3;
    case 'WM' 
%       numstim = 4;
       numstim = 8;
   otherwise
        warning('Unexpected task type.')
end


%cd(strcat('/dcs07/smart/data/HRF/bogdan_hrf/stats/fsl_noar_24vec_csf_spikes/results/', Task,'/'))
cd(strcat('/dcs07/smart/data/HRF/bogdan_hrf/stats/', Type,'_24vec_csf_spikes/results/', Task,'/'))

% Jamie_HCP_Results/stats/spm_24vec_csf_spikes_noar_spline4s/results/118528/WM/Residuals_4D_sub-118528_task-WM_dir-LR.nii.gz



names = filenames('*');

nsub = length(names);

Results = zeros(489, 41, nsub, numstim);
Flags = zeros(numstim, nsub);


atlas_obj = load_atlas('canlab2018');
%obj_template = fmri_data(strcat(names{1},'/',Direction,'/res4d.nii'));

%savename = strcat('/dcs04/smart/data/HRF_Project/Results/HRF_Resid_',Task,'_',Direction,'.mat');
savename = strcat('/dcs04/smart/data/HRF_Project/Results/', Type,'_HRF_Resid_', Task,'block_', Direction,'.mat');

for sub=1:nsub

    if ((exist(strcat(names{sub},'/',Direction,'/res4d.nii.gz'), 'file') > 0) | (exist(strcat(names{sub},'/',Direction,'/res4d.nii'), 'file') > 0))
        [Runc, Flag] = GenerateOnsetsHCP(names{sub}, Task, Direction); 
        obj = fmri_data(strcat(names{sub},'/',Direction,'/res4d.nii'));
        TC = extract_data(atlas_obj,obj);
    
        [HRF] = Fit_sFIR_all(TC', TR, Runc, T, 1, []);
        Flags(:,sub) = Flag;
    
        tmp = 1;
        for i=1:numstim
            if (Flag(i) == 1), Results(:,:,sub, i) = squeeze(HRF(:,:,tmp)); tmp=tmp+1;  end
        end
    end

    disp(sub)
    save(savename, 'Results', 'Flags')
  
end

% The result is a single 41-point-long vector that represents the
%  average shape of the residual signal following the stimulus.





