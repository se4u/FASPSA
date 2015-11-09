%{
Filename    : analyze_runs.m
Description : A script to analyze the results of `compare_algorithms.m`
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 23:29:48 2015 (-0500)
Last-Updated: .
By: .
Update #: 0

`compare_algorithms.m` produces a mat file called `sso_project.mat` that contains
the struct `results_struct` that contains keys of the form
[algorithm]_[dimension]_[runidx] for example:

Adaptive2SPSA_60_49_time_taken
Adaptive2SPSA_60_49_final_loss
Adaptive2SPSA_60_49_final_sqdist
Adaptive2SPSA_60_49_iteration_count

algorithm may be one of `Adaptive2SPSA`, `FeedbackAdaptive2SPSA`,
  `EfficientAdaptive2SPSA`, `EfficientFeedbackAdaptive2SPSA`.
dimension usually ranges between 10 to 100. and the
run ranges between 1 to 50.

This script can produce many plots depending on the cell ran.
Please see the documentation for each cell for details.
%}
%% 1. a scatter plot of the time_taken versus final_loss
% incurred with a small circle around the point that indicates the variance
% in the readings. For each dimensionality a separate scatter plot is produced.
% Each scatter plot contains 4 points corresponding to the four algorithms.
load('../res/sso_project.mat');
p = 100;
runs = 50;
budget = 2000;
algorithms = {'Adaptive2SPSA', 'FeedbackAdaptive2SPSA','EfficientAdaptive2SPSA', ...
              'EfficientFeedbackAdaptive2SPSA'};
% O-red, x-blue, square-green, diamond-black.
markups = {'or', 'xb', 'sg', 'dk'};
for algorithm_idx=1:length(algorithms)
    algorithm = algorithms{algorithm_idx};
    markup = markups{algorithm_idx};
    time_readings = nan(1, runs);
    loss_readings = nan(1, runs);
    for run_idx=1:runs
        prefix = concat_all(algorithm, p, run_idx, budget);
        time_readings(run_idx) = results_struct.([prefix, '_time_taken']);
        loss_readings(run_idx) = min(1, results_struct.([prefix, '_final_loss']));
    end
    scatter(time_readings, loss_readings, markup);
    hold on;
end
line('YData', [1 1], 'LineStyle', '-', 'LineWidth', 2, 'Color','m');
title(['Final Loss vs. Time taken for fixed budget at dimension=', ...
       num2str(p)]);
xlabel('Time');
ylabel('Loss');
ylim([0 1.1]);
legend(algorithms);
saveas(gcf, 'analyze_runs', 'png');
exit;
