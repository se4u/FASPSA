%{
Filename    : compare_algorithms.m
Description : Code to compare the Feedback|AdaptiveSPSA algorithms with their efficient counterparts.
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 17:29:17 2015 (-0500)
Last-Updated: .
By: .
Update #: 0
%}
if maxNumCompThreads() > 1
    disp('Start matlab with -singleCompThread flag.');
    exit(1);
else
    disp(computer());
end
close all; clear; clc;
rand('seed',31415927);
randn('seed',3111113);
sigma = 0.05;
sequence_param_cell.alpha = 0.602;
sequence_param_cell.gamma = 0.101;
sequence_param_cell.a_numerator = 1e-2;
sequence_param_cell.c_numerator = 1e-1;
name_fn_struct = struct();
name_fn_struct.Adaptive2SPSA = @Adaptive2SPSA;
name_fn_struct.FeedbackAdaptive2SPSA = @FeedbackAdaptive2SPSA;
name_fn_struct.EfficientAdaptive2SPSA = @EfficientAdaptive2SPSA;
name_fn_struct.EfficientFeedbackAdaptive2SPSA = ...
    @EfficientFeedbackAdaptive2SPSA;
name_fn_cell = fieldnames(name_fn_struct);
results_struct = struct();
% for multiple budgets.
% for multiple dimensions.
% for multiple runs.
% for different algorithms
% Run algorithm with dimensions.
% collate result.
% The total memory of the struct would not exceed 60MB.
for budget=[2000 10000];
for p=[10 20 30 40 50 60] % for multiple dimensions.
    init_theta = 0.2 * ones(p, 1);
    true_loss_fn = quartic_loss_factory(p);
    target_fn = noisy_function_factory(true_loss_fn, sigma);
    true_optimal_theta = zeros(p, 1);
    for run_idx=1:50 % for multiple runs.
        for name_fn_idx=1:length(name_fn_cell) % for different algorithms
            name_fn = name_fn_cell{name_fn_idx};
            FN = getfield(name_fn_struct, name_fn);
            common_prefix = concat_all(name_fn, p, run_idx, budget);
            disp(common_prefix);
            % Run the algorithm.
            [iteration_count, theta, time_taken, loss_sequence, mad_sequence] = FN(budget, target_fn, init_theta, true_loss_fn, ...
                       true_optimal_theta, sequence_param_cell);
            % collate result.
            results_struct.([common_prefix, '_time_taken']) = time_taken;
            results_struct.([common_prefix, '_final_loss']) = ...
                loss_sequence(length(loss_sequence));
            results_struct.([common_prefix, '_final_mad']) = ...
                mad_sequence(length(mad_sequence));
            results_struct.([common_prefix, '_iteration_count']) = ...
                iteration_count;
            results_struct.([common_prefix, '_loss_sequence']) = ...
                loss_sequence;
            results_struct.([common_prefix, '_mad_sequence']) = ...
                mad_sequence;
        end
    end
    % Write/Overwrite intermediate result files that can be observed separately.
    save('sso_project_intermediate.mat', 'results_struct');
end
end
save('sso_project.mat', 'results_struct');
exit;
