dat = readmatrix('data\TAFC_mean.csv');
TR = 0.72; 
load('data\WMblock_parallel_results.mat', 'TAFC');   

mean_accuracy = dat(:, 2:end);
contrast_names = {
    'body vs faces', 
    'body vs places', 
    'body vs tools', 
    'faces vs places', 
    'faces vs tools', 
    'places vs tools'
};

time_points = 1:size(mean_accuracy, 1);

%% Subplots for each contrast 
figure;
sgtitle('Forced Choice Accuracy', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:6
    subplot(3, 2, i);
    plot(time_points, mean_accuracy(:, i), 'LineWidth', 2);
   
    hold on;
    yline(0.5, 'r--', 'Chance');
    hold off;
    
    ylim([0.4, 1.0]); 
    title(contrast_names{i});
    xlabel('Time Point');
    ylabel('Accuracy');
    grid on;
end

%% Overall Discriminability Ranking
figure;
overall_mean = mean(mean_accuracy, 1);
[sorted_means, sort_index] = sort(overall_mean, 'descend');
sorted_names = contrast_names(sort_index);
b = bar(sorted_means);
set(gca, 'XTickLabel', sorted_names);
xtickangle(45); 
hold on;
yline(0.5, 'r--', 'Chance Level');
hold off;

title('Overall Average Forced Choice Accuracy by Contrast', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Average Accuracy (across all time points)');
ylim([0.4, 1.0]);

text(1:6, sorted_means, num2str(sorted_means', '%.2f'), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');  
grid on;


% Display the ranked list
fprintf('%-4s | %-20s | %s\n', 'Rank', 'Contrast', 'Mean Accuracy');
disp(repmat('-', 1, 45));
for i = 1:length(sorted_names)
    fprintf('#%-3d | %-20s |    %.4f\n', ...
        i, sorted_names{i}, sorted_means(i));
end
fprintf('\n');

%% Correlation matrix
R = corr(mean_accuracy);
disp('Trajectory Correlation Matrix (R):');
disp(R);

figure;
h = heatmap(contrast_names, contrast_names, R);
h.Title = 'Similarity (Correlation) of Trajectories';
h.XLabel = 'Contrast';
h.YLabel = 'Contrast';

h.ColorLimits = [-1, 1]; 


%%  Peak Accuracy 
[peak_acc, peak_time_index] = max(mean_accuracy, [], 1);
peak_time_sec = (peak_time_index - 1) * TR;

fprintf('%-20s | %s | %s\n', 'Contrast', 'Peak Accuracy', 'Time-to-Peak (s)');
disp(repmat('-', 1, 60));
for i = 1:length(contrast_names)
    fprintf('%-20s |      %.3f      |     %.2f s\n', ...
        contrast_names{i}, peak_acc(i), peak_time_sec(i));
end
fprintf('\n');

%% --- 4. Statistical Comparison of Separability ---
mean_sep_bootstraps = squeeze(mean(TAFC, 1));

% --- STEP 3: Run the paired t-test ---
% Example: 'faces vs places' (contrast 4) vs 'body vs tools' (contrast 3)
contrast_A_data = mean_sep_bootstraps(:, 4); % Faces vs Places
contrast_B_data = mean_sep_bootstraps(:, 3); % Body vs Tools

[h, p, ci, stats] = ttest(contrast_A_data, contrast_B_data);

fprintf('Statistical Test: "Faces vs Places" > "Body vs Tools"\n');
fprintf('p-value: %.5f\n', p);
fprintf('t-stat:  %.3f\n', stats.tstat);

if h == 1 && stats.tstat > 0
     disp('Result: "Faces vs Places" is significantly MORE separable (p < 0.05).');
elseif h == 1 && stats.tstat < 0
     disp('Result: "Body vs Tools" is significantly MORE separable (p < 0.05).');
else
     disp('Result: The difference is NOT statistically significant.');
end

