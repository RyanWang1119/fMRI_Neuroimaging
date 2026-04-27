% For the specified name, task, and direction. Extract the onset file and 
% convert the timings to the "fMRI-data" domain of TRs. Return a cell array

function [Runc, Flag] = GenerateOnsetsHCP(name, task, direction)
%
% Generate onsets for HCP tasks

TR = 0.72;
path = strcat('/dcs07/smart/data/human-connectome-project-openaccess/HCP1200/',name,'/MNINonLinear/Results/');


%        cd '/Users/martinlindquist/Desktop/Data/dcs04/legacy-dcs01-oasis/hpc/disk1/113619/MNINonLinear/Results'
%        cd(strcat('tfMRI_',task,'_',direction,'/EVs'))

switch task

    case 'EMOTION' 

        % Generate onsets for HCP Emotion task 
        nvol = 176;         
        path_emotion = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');

        S1 = ExtractCondition(strcat(path_emotion,'fear.txt'), TR, nvol);
        S2 = ExtractCondition(strcat(path_emotion,'neut.txt'), TR, nvol);
        
        Runc{1} = S1;
        Runc{2} = S2;

        tmp = 1;
        Flag = zeros(2,1);
        if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
        if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; end

    case 'GAMBLING' 

        % Generate onsets for HCP GAMBLING task
     	nvol = 253;
        path_gambling = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');

        S1 = ExtractCondition(strcat(path_gambling,'loss_event.txt'), TR, nvol);
        S2 = ExtractCondition(strcat(path_gambling,'neut_event.txt'), TR, nvol);
        S3 = ExtractCondition(strcat(path_gambling,'win_event.txt'), TR, nvol);
        
        tmp = 1;
        Flag = zeros(3,1);
        if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
        if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; tmp = tmp+1; end
        if (S3 ~= -1), Runc{tmp} = S3(1:nvol); Flag(3) =1; end   

    case 'LANGUAGE' 

        % Generate onsets for HCP LANGUAGE task
        nvol = 316;
        path_language = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');

        S1 = ExtractCondition(strcat(path_language,'present_math.txt'), TR, nvol);
        S2 = ExtractCondition(strcat(path_language,'question_math.txt'), TR, nvol);
        S3 = ExtractCondition(strcat(path_language,'response_math.txt'), TR, nvol);
        S4 = ExtractCondition(strcat(path_language,'present_story.txt'), TR, nvol);
        S5 = ExtractCondition(strcat(path_language,'question_story.txt'), TR, nvol);
        S6 = ExtractCondition(strcat(path_language,'response_story.txt'), TR, nvol);
 
        tmp = 1;
        Flag = zeros(6,1);
        if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
        if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; tmp = tmp+1; end
        if (S3 ~= -1), Runc{tmp} = S3(1:nvol); Flag(3) =1; tmp = tmp+1; end
        if (S4 ~= -1), Runc{tmp} = S4(1:nvol); Flag(4) =1; tmp = tmp+1; end
        if (S5 ~= -1), Runc{tmp} = S5(1:nvol); Flag(5) =1; tmp = tmp+1; end
        if (S6 ~= -1), Runc{tmp} = S6(1:nvol); Flag(6) =1; end

   case 'MOTOR'

        % Generate onsets for HCP MOTOR task
        
%        cd '/Users/martinlindquist/Desktop/Data/dcs04/legacy-dcs01-oasis/hpc/disk1/113619/MNINonLinear/Results'
%        cd(strcat('tfMRI_',task,'_',direction,'/EVs'))


%        cd '/Users/martinlindquist/Library/CloudStorage/Dropbox/MOTOR_DATA/125525/tfMRI_MOTOR_RL/EVs'
        nvol = 284;
        path_motor = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');

        S1 = ExtractCondition(strcat(path_motor,'cue.txt'), TR, nvol);
        S2 = ExtractCondition(strcat(path_motor,'lf.txt'), TR, nvol);
        S3 = ExtractCondition(strcat(path_motor,'rf.txt'), TR, nvol);
        S4 = ExtractCondition(strcat(path_motor,'lh.txt'), TR, nvol);
        S5 = ExtractCondition(strcat(path_motor,'rh.txt'), TR, nvol);
        S6 = ExtractCondition(strcat(path_motor,'t.txt'), TR, nvol);
                
        tmp = 1;
        Flag = zeros(6,1);
        if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
        if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; tmp = tmp+1; end
        if (S3 ~= -1), Runc{tmp} = S3(1:nvol); Flag(3) =1; tmp = tmp+1; end
        if (S4 ~= -1), Runc{tmp} = S4(1:nvol); Flag(4) =1; tmp = tmp+1; end
        if (S5 ~= -1), Runc{tmp} = S5(1:nvol); Flag(5) =1; tmp = tmp+1; end
        if (S6 ~= -1), Runc{tmp} = S6(1:nvol); Flag(6) =1; end

    case 'RELATIONAL' 
 
        % Generate onsets for HCP Relational task
          
%        cd(strcat('/Users/martinlindquist/Library/CloudStorage/Dropbox/WM_DATA/',name,'/tfMRI_WM_RL/EVs'))
        
        nvol = 232;
        path_rel = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');
       
        S1 = ExtractCondition(strcat(path_rel,'error.txt'), TR, nvol);
        S2 = ExtractCondition(strcat(path_rel,'match.txt'), TR, nvol);
        S3 = ExtractCondition(strcat(path_rel,'relation.txt'), TR, nvol);

        tmp = 1;
        Flag = zeros(3,1);
        if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
        if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; tmp = tmp+1; end
        if (S3 ~= -1), Runc{tmp} = S3(1:nvol); Flag(3) =1; end
% 
%     case 'WM' 
% 
%         % Generate onsets for HCP WM task
% 
% %        cd(strcat('/Users/martinlindquist/Library/CloudStorage/Dropbox/WM_DATA/',name,'/tfMRI_WM_RL/EVs'))
% 	    nvol = 405;
%         path_wm = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');
% 
%         S1 = ExtractCondition(strcat(path_wm,'0bk_cor.txt'), TR, nvol);
%         S2 = ExtractCondition(strcat(path_wm,'0bk_err.txt'), TR, nvol);
%         S3 = ExtractCondition(strcat(path_wm,'2bk_cor.txt'), TR, nvol);
%         S4 = ExtractCondition(strcat(path_wm,'2bk_err.txt'), TR, nvol);
% 
%         tmp = 1;
%         Flag = zeros(4,1);
%         if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
%         if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; tmp = tmp+1; end
%         if (S3 ~= -1), Runc{tmp} = S3(1:nvol); Flag(3) =1; tmp = tmp+1; end
%         if (S4 ~= -1), Runc{tmp} = S4(1:nvol); Flag(4) =1; end


    case 'WM' 

        % Generate onsets for HCP WM task
      
%        cd(strcat('/Users/martinlindquist/Library/CloudStorage/Dropbox/WM_DATA/',name,'/tfMRI_WM_RL/EVs'))
	    nvol = 405;
        path_wm = strcat(path,'tfMRI_',task,'_',direction,'/EVs/');
    
        S1 = ExtractCondition(strcat(path_wm,'0bk_body.txt'), TR, nvol);
        S2 = ExtractCondition(strcat(path_wm,'0bk_faces.txt'), TR, nvol);
        S3 = ExtractCondition(strcat(path_wm,'0bk_places.txt'), TR, nvol);
        S4 = ExtractCondition(strcat(path_wm,'0bk_tools.txt'), TR, nvol);
        S5 = ExtractCondition(strcat(path_wm,'2bk_body.txt'), TR, nvol);
        S6 = ExtractCondition(strcat(path_wm,'2bk_faces.txt'), TR, nvol);
        S7 = ExtractCondition(strcat(path_wm,'2bk_places.txt'), TR, nvol);
        S8 = ExtractCondition(strcat(path_wm,'2bk_tools.txt'), TR, nvol);
         
        tmp = 1;
        Flag = zeros(4,1);
        if (S1 ~= -1), Runc{tmp} = S1(1:nvol); Flag(1) =1; tmp = tmp+1; end
        if (S2 ~= -1), Runc{tmp} = S2(1:nvol); Flag(2) =1; tmp = tmp+1; end
        if (S3 ~= -1), Runc{tmp} = S3(1:nvol); Flag(3) =1; tmp = tmp+1; end
        if (S4 ~= -1), Runc{tmp} = S4(1:nvol); Flag(4) =1; tmp = tmp+1; end
        if (S5 ~= -1), Runc{tmp} = S5(1:nvol); Flag(5) =1; tmp = tmp+1; end
        if (S6 ~= -1), Runc{tmp} = S6(1:nvol); Flag(6) =1; tmp = tmp+1; end
        if (S7 ~= -1), Runc{tmp} = S7(1:nvol); Flag(7) =1; tmp = tmp+1; end
        if (S8 ~= -1), Runc{tmp} = S8(1:nvol); Flag(8) =1; end

   otherwise
        warning('Unexpected task type.')
end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sub Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function S = ExtractCondition(cond, TR, nvol)

cnd = load(cond);
S = -1;
if (~isempty(cnd))
    cnd = round(cnd/TR);
    cnd(:,1) = max(cnd(:,1),1);
    num = size(cnd,1);
    S = zeros(nvol,1);
    for i=1:num
        S(cnd(i,1):(cnd(i,1) + cnd(i,2))) = 1;
    end
end



end

