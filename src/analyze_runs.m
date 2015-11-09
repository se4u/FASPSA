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
clear; clc; close all;
load('../res/sso_project.mat');
p = 10;
runs = 2;
budget = 25000;
algorithms = {'Adaptive2SPSA', 'FeedbackAdaptive2SPSA','EfficientAdaptive2SPSA', ...
              'EfficientFeedbackAdaptive2SPSA'};
% O-red, x-blue, square-green, diamond-black.
markups = {'or', 'xb', 'sg', 'dk'};
for algorithm_idx=1:length(algorithms)
    algorithm = algorithms{algorithm_idx}
    markup = markups{algorithm_idx};
    time_readings = nan(1, runs);
    loss_readings = nan(1, runs);
    errtheta_readings = nan(1, runs);
    losssq_readings = nan(1, runs);
    errtheta_init = nan;
    losssq_init = nan;
    for run_idx=1:runs
        prefix = concat_all(algorithm, p, run_idx, budget);
        time_taken = results_struct.([prefix, '_time_taken']);
        loss_seq = results_struct.([prefix, '_loss_sequence']);
        sqdist_seq = results_struct.([prefix, '_sqdist_sequence']);

        time_readings(run_idx) = time_taken;
        loss_readings(run_idx) = loss_seq(length(loss_seq));
        errtheta_readings(run_idx) = sqdist_seq(length(sqdist_seq));
        losssq_readings(run_idx) = loss_seq(length(loss_seq))^2;

        if isnan(losssq_init)
            losssq_init = loss_seq(1);
        else
            assert(losssq_init == loss_seq(1));
        end
        if isnan(errtheta_init)
            errtheta_init = sqdist_seq(1);
        else
            assert(errtheta_init == sqdist_seq(1));
        end
    end
    norm_sqloss = mean(losssq_readings)/losssq_init
    norm_theta = sqrt(mean(errtheta_readings)/errtheta_init)
    scatter(time_readings, loss_readings, markup);
    hold on;
end
title(['Final Loss vs. Time taken for fixed budget at dimension=', ...
       num2str(p)]);
xlabel('Time');
ylabel('Loss');
legend(algorithms, 'Location','NorthEastOutside');
saveas(gcf, '../res/analyze_runs', 'png');
