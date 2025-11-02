load('data\HRF_Resid_WMblock_LR.mat')

time_points = 41;
regions = 489;
subjects = 393;

Cond = [[5 1 6 2]; [5 1 7 3]; [5 1 8 4]; [6 2 7 3]; [6 2 8 4]; [7 3 8 4]];
% [5 1 6 2]: (2bk_body - 0bk_body) vs. (2bk_faces - 0bk_faces)
% [5 1 7 3]: (2bk_body - 0bk_body) vs. (2bk_places - 0bk_places)
% [5 1 8 4]: (2bk_body - 0bk_body) vs. (2bk_tools - 0bk_tools)
% [6 2 7 3]: (2bk_faces - 0bk_faces) vs. (2bk_places - 0bk_places)
% [6 2 8 4]: (2bk_faces - 0bk_faces) vs. (2bk_tools - 0bk_tools)
% [7 3 8 4]: (2bk_places - 0bk_places) vs. (2bk_tools - 0bk_tools)

bsamp = 20;

numWorkers = 8;  
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

TA = zeros(time_points, bsamp, 6);
TAFC = zeros(time_points, bsamp, 6);
Beta = zeros(time_points, bsamp, 6, regions);

parfor rep = 1:bsamp
    fprintf('Starting bootstrap sample %d/%d\n', rep, bsamp);
    
    TA_local = zeros(time_points, 6);
    TAFC_local = zeros(time_points, 6);
    Beta_local = zeros(time_points, 6, regions);

    for contrast = 1:6
        for t = 1:time_points
            
            X1 = (squeeze(Results(:,t,:,Cond(contrast,1))) - squeeze(Results(:,t,:,Cond(contrast,2))))';
            X2 = (squeeze(Results(:,t,:,Cond(contrast,3))) - squeeze(Results(:,t,:,Cond(contrast,4))))';
            
            Y = zeros(2*subjects, 1);
            Y(1:subjects) = 1;
            Y(subjects+1:2*subjects) = -1;
            ids = [(1:subjects)'; (1:subjects)'];
            
            train = [1:314 (subjects+1):(subjects+314)]';
            test = [315:subjects (subjects+315):(2*subjects)]';
            
            % Bootstrap sampling
            bs_sample = randsample(314, 314, true);
            Xtrain_bs = [X1(bs_sample,:); X2(bs_sample,:)];
            Xtest = [X1(test(1:length(test)/2),:); X2(test(1:length(test)/2),:)];
            
            Ytrain = Y(train);
            Ytest = Y(test);
            idstrain = ids(train);
            
            S = xval_SVM(Xtrain_bs, Ytrain, idstrain, 'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot');
            
            beta_weights = S.ClassificationModel.Beta; 
            Beta_local(t, contrast, :) = beta_weights';

            [yhat, score] = predict(S.ClassificationModel, Xtest);
            TA_local(t, contrast) = sum(Ytest == yhat) / length(Ytest);
           
            pairedscores = reshape(score(:, 1), length(Ytest)/2, 2);
            iscorrect = diff(pairedscores') > 0;
            TAFC_local(t, contrast) = sum(iscorrect) / length(iscorrect);

        end
    end
    
    TA(:, rep, :) = TA_local;
    TAFC(:, rep, :) = TAFC_local;
    Beta(:, rep, :, :) = Beta_local;
    fprintf('Completed bootstrap sample %d/%d\n', rep, bsamp);
end

TA_mean = squeeze(mean(TA, 2));
TA_ste = squeeze(std(TA, 0, 2)) / sqrt(bsamp);
TAFC_mean = squeeze(mean(TAFC, 2));
TAFC_ste = squeeze(std(TAFC, 0, 2)) / sqrt(bsamp);
Beta_mean = squeeze(mean(Beta, 2));
Beta_std = squeeze(std(Beta, 0, 2)); 

%% Plotting
figure
subplot(2, 1, 1)
plot(TA_mean, 'LineWidth', 2)
ylim([0 1])
title('Test Accuracy')
legend('body vs faces', 'body vs places', 'body vs tools', 'faces vs places', 'faces vs tools', 'places vs tools')
xlabel('Time Point')
ylabel('Accuracy')

subplot(2, 1, 2)
plot(TAFC_mean, 'LineWidth', 2)
ylim([0 1])
title('Forced Choice Accuracy')
legend('body vs faces', 'body vs places', 'body vs tools', 'faces vs places', 'faces vs tools', 'places vs tools')
xlabel('Time Point')
ylabel('Accuracy')

% Save results
save('WMblock_parallel_results.mat', 'TA', 'TAFC', 'TA_mean', 'TA_ste', ...
     'TAFC_mean', 'TAFC_ste', 'Beta', 'Beta_mean', 'Beta_std');
delete(gcp('nocreate'));

%% save
contrast_names = {'body_vs_faces', 'body_vs_places', 'body_vs_tools', ...
                  'faces_vs_places', 'faces_vs_tools', 'places_vs_tools'};

% Test Accuracy Mean
TA_mean_table = array2table(TA_mean, 'VariableNames', contrast_names);
TA_mean_table.TimePoint = (1:time_points)';
TA_mean_table = movevars(TA_mean_table, 'TimePoint', 'Before', 1);
writetable(TA_mean_table, 'TA_mean.csv');

% Test Accuracy Standard Error
TA_ste_table = array2table(TA_ste, 'VariableNames', contrast_names);
TA_ste_table.TimePoint = (1:time_points)';
TA_ste_table = movevars(TA_ste_table, 'TimePoint', 'Before', 1);
writetable(TA_ste_table, 'TA_ste.csv');

% Forced Choice Accuracy Mean
TAFC_mean_table = array2table(TAFC_mean, 'VariableNames', contrast_names);
TAFC_mean_table.TimePoint = (1:time_points)';
TAFC_mean_table = movevars(TAFC_mean_table, 'TimePoint', 'Before', 1);
writetable(TAFC_mean_table, 'TAFC_mean.csv');

% Forced Choice Accuracy Standard Error
TAFC_ste_table = array2table(TAFC_ste, 'VariableNames', contrast_names);
TAFC_ste_table.TimePoint = (1:time_points)';
TAFC_ste_table = movevars(TAFC_ste_table, 'TimePoint', 'Before', 1);
writetable(TAFC_ste_table, 'TAFC_ste.csv');