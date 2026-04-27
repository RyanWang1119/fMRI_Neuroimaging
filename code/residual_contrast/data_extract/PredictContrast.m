load('HRF_Resid_WMblock_LR.mat')


time_points = 41;
regions = 489;
subjects = 393;


Y = zeros(2*subjects,1);
X = zeros(2*subjects, regions);

accuracy = zeros(time_points,6);
precision = zeros(time_points,6);
recall = zeros(time_points,6);
f1_score = zeros(time_points,6);

test_acc = zeros(time_points,6);
test_acc_forcedchoice = zeros(time_points,6);

Cond = [[5 1 6 2]; [5 1 7 3]; [5 1 8 4]; [6 2 7 3]; [6 2 8 4]; [7 3 8 4]];

% [5 1 6 2]: (2bk_body - 0bk_body) vs. (2bk_faces - 0bk_faces)
% [5 1 7 3]: (2bk_body - 0bk_body) vs. (2bk_places - 0bk_places)
% [5 1 8 4]: (2bk_body - 0bk_body) vs. (2bk_tools - 0bk_tools)
% [6 2 7 3]: (2bk_faces - 0bk_faces) vs. (2bk_places - 0bk_places)
% [6 2 8 4]: (2bk_faces - 0bk_faces) vs. (2bk_tools - 0bk_tools)
% [7 3 8 4]: (2bk_places - 0bk_places) vs. (2bk_tools - 0bk_tools)
% 2bk - 0bk cancels out the nuisance effects and extract the WM features.
% WM*body vs WM*face.

bsamp = 20;

TA = zeros(time_points,bsamp,6);
TAFC = zeros(time_points,bsamp,6);

for rep = 1:bsamp,
for contrast=1:6
    for t = 1:time_points

        X1 = (squeeze(Results(:,t,:,Cond(contrast,1))) - squeeze(Results(:,t,:,Cond(contrast,2))))';
        X2 = (squeeze(Results(:,t,:,Cond(contrast,3))) - squeeze(Results(:,t,:,Cond(contrast,4))))';  
        X = [X1; X2];
        
        Y(1: subjects) = 1;
        Y(subjects+1: 2*subjects) = -1;
        ids = [(1:subjects)'; (1:subjects)']; 
        
        train = [1:314 (subjects+1):(subjects+314)]';
        test = [315:subjects (subjects+315):(2*subjects)]';

        Xtrain = X(train,:);


        bs_sample = randsample(314, 314, true);
        Xtrain_bs = [X1(bs_sample,:); X2(bs_sample,:)];

        Xtest = X(test,:);

        Ytrain = Y(train);
        Ytest = Y(test);

        idstrain = ids(train);
        idstest = ids(test);

        S = xval_SVM(Xtrain_bs, Ytrain, idstrain, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');
  
        cm = confusionchart(S.Y, S.yfit);
        C = cm.NormalizedValues;
        TN = C(1,1);
        FP = C(1,2);
        FN = C(2,1);
        TP = C(2,2);
        
        % Compute metrics
        accuracy(t,contrast) = (TP + TN) / sum(C(:));
        precision(t,contrast) = TP / (TP + FP);
        recall(t,contrast) = TP / (TP + FN);
        f1_score(t,contrast) = 2 * (precision(t,contrast) * recall(t,contrast)) / (precision(t,contrast) + recall(t,contrast));
    
        disp(t);

        [yhat, score] = predict(S.ClassificationModel, Xtest);
        test_acc(t,contrast) = sum(Ytest == yhat) ./ length(Ytest);
        TA(t,rep,contrast) = sum(Ytest == yhat) ./ length(Ytest);


        pairedscores = reshape(score(:, 1), length(Ytest)/2, 2);
        iscorrect = diff(pairedscores') > 0;
        test_acc_forcedchoice(t,contrast) = sum(iscorrect) ./ length(iscorrect);
        TAFC(t,rep,contrast) = sum(iscorrect) ./ length(iscorrect);


    end

end

end


TA_mean = zeros(time_points,6);
TA_ste = zeros(time_points,6);
TAFC_mean = zeros(time_points,6);
TAFC_ste = zeros(time_points,6);


for contrast=1:6
    for t = 1:time_points

        TA_mean(t,contrast) = mean(squeeze(TA(t,:,contrast)));
        TA_ste(t,contrast) = ste(squeeze(TA(t,:,contrast))');
        TAFC_mean(t,contrast) = mean(squeeze(TA(t,:,contrast)));
        TAFC_ste(t,contrast) = ste(squeeze(TA(t,:,contrast))');
       
    end
end




figure
subplot 211
plot(test_acc,LineWidth=2)
ylim([0 1])
legend('body vs faces', 'body vs places', 'body vs tools', 'faces vs places', 'faces vs tools', 'places vs tools')
subplot 212
plot(test_acc_forcedchoice,LineWidth=2)
ylim([0 1])
legend('body vs faces', 'body vs places', 'body vs tools', 'faces vs places', 'faces vs tools', 'places vs tools')

%%

% Contrasts:
% 1. 'LeftFoot_v_LeftHand' 2 v 4
% 2. 'RightFoot_v_RightHand' 3 v 5
% 3. 'LeftHand_v_RightHand' 4 v 5
% 4. 'Motor_tongue_v_cue' 6 v 1

Cond = [[2 4 3 5]; [2 4 4 5]; [2 4 6 1]; [3 5 4 5]; [3 5 6 1]; [4 5 6 1]];
    

% S1 = ExtractCondition(strcat(path_motor,'cue.txt'), TR, nvol);
% S2 = ExtractCondition(strcat(path_motor,'lf.txt'), TR, nvol);
% S3 = ExtractCondition(strcat(path_motor,'rf.txt'), TR, nvol);
% S4 = ExtractCondition(strcat(path_motor,'lh.txt'), TR, nvol);
% S5 = ExtractCondition(strcat(path_motor,'rh.txt'), TR, nvol);
% S6 = ExtractCondition(strcat(path_motor,'t.txt'), TR, nvol);

load('HRF_Resid_MOTOR_LR.mat')


time_points = 41;
regions = 489;
subjects = 405;


Y = zeros(2*subjects,1);
X = zeros(2*subjects, regions);

accuracy = zeros(time_points,6);
precision = zeros(time_points,6);
recall = zeros(time_points,6);
f1_score = zeros(time_points,6);

test_acc = zeros(time_points,6);
test_acc_forcedchoice = zeros(time_points,6);

for contrast=1:6
    for t = 1:time_points

        X1 = (squeeze(Results(:,t,:,Cond(contrast,1))) - squeeze(Results(:,t,:,Cond(contrast,2))))';
        X2 = (squeeze(Results(:,t,:,Cond(contrast,3))) - squeeze(Results(:,t,:,Cond(contrast,4))))';  
        X = [X1; X2];
        
        Y(1: subjects) = 1;
        Y(subjects+1: 2*subjects) = -1;
        ids = [(1:subjects)'; (1:subjects)']; 
        
        train = [1:314 (subjects+1):(subjects+314)]';
        test = [315:subjects (subjects+315):(2*subjects)]';

        Xtrain = X(train,:);
        Xtest = X(test,:);

        Ytrain = Y(train);
        Ytest = Y(test);

        idstrain = ids(train);
        idstest = ids(test);

        S = xval_SVM(Xtrain, Ytrain, idstrain, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');
  
        cm = confusionchart(S.Y, S.yfit);
        C = cm.NormalizedValues;
        TN = C(1,1);
        FP = C(1,2);
        FN = C(2,1);
        TP = C(2,2);
        
        % Compute metrics
        accuracy(t,contrast) = (TP + TN) / sum(C(:));
        precision(t,contrast) = TP / (TP + FP);
        recall(t,contrast) = TP / (TP + FN);
        f1_score(t,contrast) = 2 * (precision(t,contrast) * recall(t,contrast)) / (precision(t,contrast) + recall(t,contrast));
    
        disp(t);

        [yhat, score] = predict(S.SVMModel, Xtest);
        test_acc(t,contrast) = sum(Ytest == yhat) ./ length(Ytest);


        pairedscores = reshape(score(:, 1), length(Ytest)/2, 2);
        iscorrect = diff(pairedscores') > 0;
        test_acc_forcedchoice(t,contrast) = sum(iscorrect) ./ length(iscorrect);


    end

end

% 1. 'LeftFoot_v_LeftHand' 2 v 4
% 2. 'RightFoot_v_RightHand' 3 v 5
% 3. 'LeftHand_v_RightHand' 4 v 5
% 4. 'Motor_tongue_v_cue' 6 v 

figure
subplot 211
plot(test_acc,LineWidth=2)
ylim([0 1])
legend('LF-LH vs RF-RH', 'LF-LH vs LH-RH', 'LF-LH vs T-Cue', 'RF-RH vs LH-RH', 'RF-RH vs T-Cue', 'LH-RH vs T-Cue')
subplot 212
plot(test_acc_forcedchoice,LineWidth=2)
ylim([0 1])
legend('LF-LH vs RF-RH', 'LF-LH vs LH-RH', 'LF-LH vs T-Cue', 'RF-RH vs LH-RH', 'RF-RH vs T-Cue', 'LH-RH vs T-Cue')


%%


Cond = [[2 3]; [2 4]; [2 5]; [2 6]; [3 4]; [3 5]; [3 6]; [4 5]; [4 6] ; [5 6]];
    

% S1 = ExtractCondition(strcat(path_motor,'cue.txt'), TR, nvol);
% S2 = ExtractCondition(strcat(path_motor,'lf.txt'), TR, nvol);
% S3 = ExtractCondition(strcat(path_motor,'rf.txt'), TR, nvol);
% S4 = ExtractCondition(strcat(path_motor,'lh.txt'), TR, nvol);
% S5 = ExtractCondition(strcat(path_motor,'rh.txt'), TR, nvol);
% S6 = ExtractCondition(strcat(path_motor,'t.txt'), TR, nvol);

load('HRF_Resid_MOTOR_LR.mat')


time_points = 41;
regions = 489;
subjects = 405;


Y = zeros(2*subjects,1);
X = zeros(2*subjects, regions);

accuracy = zeros(time_points,10);
precision = zeros(time_points,10);
recall = zeros(time_points,10);
f1_score = zeros(time_points,10);

test_acc = zeros(time_points,10);
test_acc_forcedchoice = zeros(time_points,10);

for contrast=1:10
    for t = 1:time_points

        X1 = (squeeze(Results(:,t,:,Cond(contrast,1))))';
        X2 = (squeeze(Results(:,t,:,Cond(contrast,2))))';  
        X = [X1; X2];
        
        Y(1: subjects) = 1;
        Y(subjects+1: 2*subjects) = -1;
        ids = [(1:subjects)'; (1:subjects)']; 
        
        train = [1:314 (subjects+1):(subjects+314)]';
        test = [315:subjects (subjects+315):(2*subjects)]';

        Xtrain = X(train,:);
        Xtest = X(test,:);

        Ytrain = Y(train);
        Ytest = Y(test);

        idstrain = ids(train);
        idstest = ids(test);

        S = xval_SVM(Xtrain, Ytrain, idstrain, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');
  
        cm = confusionchart(S.Y, S.yfit);
        C = cm.NormalizedValues;
        TN = C(1,1);
        FP = C(1,2);
        FN = C(2,1);
        TP = C(2,2);
        
        % Compute metrics
        accuracy(t,contrast) = (TP + TN) / sum(C(:));
        precision(t,contrast) = TP / (TP + FP);
        recall(t,contrast) = TP / (TP + FN);
        f1_score(t,contrast) = 2 * (precision(t,contrast) * recall(t,contrast)) / (precision(t,contrast) + recall(t,contrast));
    
        disp(t);

        [yhat, score] = predict(S.SVMModel, Xtest);
        test_acc(t,contrast) = sum(Ytest == yhat) ./ length(Ytest);


        pairedscores = reshape(score(:, 1), length(Ytest)/2, 2);
        iscorrect = diff(pairedscores') > 0;
        test_acc_forcedchoice(t,contrast) = sum(iscorrect) ./ length(iscorrect);


    end

end

% S2 = ExtractCondition(strcat(path_motor,'lf.txt'), TR, nvol);
% S3 = ExtractCondition(strcat(path_motor,'rf.txt'), TR, nvol);
% S4 = ExtractCondition(strcat(path_motor,'lh.txt'), TR, nvol);
% S5 = ExtractCondition(strcat(path_motor,'rh.txt'), TR, nvol);
% S6 = ExtractCondition(strcat(path_motor,'t.txt'), TR, nvol);

figure
subplot 211
plot(test_acc,LineWidth=2)
ylim([0 1])
legend('LF vs RF', 'LF vs LH', 'LF vs RH', 'LF vs T', 'RF vs LH', 'RF vs RH', 'RF Vs T', 'LH vs RH', 'LH Vs T', 'RH vs T');  
subplot 212
plot(test_acc_forcedchoice,LineWidth=2)
ylim([0 1])
legend('LF vs RF', 'LF vs LH', 'LF vs RH', 'LF vs T', 'RF vs LH', 'RF vs RH', 'RF Vs T', 'LH vs RH', 'LH Vs T', 'RH vs T');  





% %%
% 
% 
% load('HRF_Resid_WMblock_LR.mat')
% 
% nrep = 3;
% 
% accuracy = zeros(time_points,nrep);
% precision = zeros(time_points,nrep);
% recall = zeros(time_points,nrep);
% f1_score = zeros(time_points,nrep);
% 
% test_acc = zeros(time_points,nrep);
% test_acc_forcedchoice = zeros(time_points,nrep);
% 
% 
% for i = 1:nrep
% 
%     contrast=1;
% 
%     Y(1: subjects) = 1;
%     Y(subjects+1: 2*subjects) = -1;
%     ids = [(1:subjects)'; (1:subjects)']; 
% 
%     train = [1:314 (subjects+1):(subjects+314)]';
%     test = [315:subjects (subjects+315):(2*subjects)]';
%  %   train_perm = randperm(314);
% 
% 
%     Ytrain = Y(train);
%     Ytest = Y(test);
% 
%     idstrain = ids(train);
%     idstest = ids(test);
% 
%     for t = 1:time_points
% 
%         X1 = (squeeze(Results(:,t,:,Cond(contrast,1))) - squeeze(Results(:,t,:,Cond(contrast,2))))';
%         X2 = (squeeze(Results(:,t,:,Cond(contrast,3))) - squeeze(Results(:,t,:,Cond(contrast,4))))';
% 
%         X = [X1; X2];
% %        Xtrain  = X(train,:);
%         Xtest = X(test,:);
%         Xtrain_perm = [X1(train_perm,:); X2(train_perm,:)];
%         idstrain_perm = [train_perm'; train_perm'];
% 
%         S = xval_SVM(Xtrain_perm, Ytrain, idstrain_perm, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');
% 
% %        S = xval_SVM(Xtrain, Ytrain, idstrain, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');
% 
%         cm = confusionchart(S.Y, S.yfit);
%         C = cm.NormalizedValues;
%         TN = C(1,1);
%         FP = C(1,2);
%         FN = C(2,1);
%         TP = C(2,2);
% 
%         % Compute metrics
%         accuracy(t,i) = (TP + TN) / sum(C(:));
%         precision(t,i) = TP / (TP + FP);
%         recall(t,i) = TP / (TP + FN);
%         f1_score(t,i) = 2 * (precision(t,i) * recall(t,i)) / (precision(t,i) + recall(t,i));
% 
%         disp(t);
% 
%         [yhat, score] = predict(S.SVMModel, Xtest);
%         test_acc(t,i) = sum(Ytest == yhat) ./ length(Ytest);
% 
%         pairedscores = reshape(score(:, 1), length(Ytest)/2, 2);
%         iscorrect = diff(pairedscores') > 0;
%         test_acc_forcedchoice(t,i) = sum(iscorrect) ./ length(iscorrect);
% 
%     end
% end
