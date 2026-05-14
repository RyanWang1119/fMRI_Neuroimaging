function plot_bootstrap_diff_curve(curve_out, timeVec, ttl)
% PLOT_BOOTSTRAP_DIFF_CURVE
%
% curve_out should be output from paired_bootstrap_curve()

    if nargin < 2 || isempty(timeVec)
        timeVec = 1:numel(curve_out.diff_mean);
    end
    if nargin < 3
        ttl = 'Bootstrap difference curve';
    end

    dm = curve_out.diff_mean(:);
    pw = curve_out.pointwise_ci;
    sb = curve_out.simul_band;

    figure; hold on;
    fill([timeVec(:); flipud(timeVec(:))], ...
         [pw(:,1); flipud(pw(:,2))], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none');
    plot(timeVec, dm, 'k', 'LineWidth', 2);
    plot(timeVec, sb(:,1), '--', 'LineWidth', 1);
    plot(timeVec, sb(:,2), '--', 'LineWidth', 1);
    yline(0, ':');
    xlabel('Time');
    ylabel('Difference');
    title(ttl);
    legend({'Pointwise 95% CI', 'Mean difference', 'Simultaneous band'}, ...
           'Location', 'best');
    box on;
end