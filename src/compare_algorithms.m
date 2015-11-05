%{
Filename    : compare_algorithms.m
Description : Code to compare the Feedback|AdaptiveSPSA algorithms with their efficient counterparts.
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 17:29:17 2015 (-0500)
Last-Updated: .
By: .
Update #: 0
%}
close all; clear; clc;
warning('off','MATLAB:maxNumCompThreads:Deprecated');
if maxNumCompThreads() > 1
    disp('Start matlab with -singleCompThread flag.');
    exit(1);
end
dir_prefix = '../res';
if exist(dir_prefix, 'dir') ~= 7
    disp(['Couldnt find directory ' dir_prefix]);
    dir_prefix = '';
    disp(['Saving files to current directory' pwd()]);
end
rand('seed',31415927);
randn('seed',3111113);
sigma = 0.05;
sequence_param_struct.alpha = 0.602;
sequence_param_struct.gamma = 0.101;
% a_numerator should be tuned. In section VIII of the 2009 paper
% professor spall said that he used a = 100.
sequence_param_struct.a_numerator = 100;
% set c to be equal to the std of the noise.
sequence_param_struct.c_numerator = sigma;
% Set c_tilda to be slightly higher than c_tilda.
sequence_param_struct.c_tilda_k_multiplier = 1.1;
sequence_param_struct.use_greedy_algorithm_b = 1;
sequence_param_struct.greedy_algorithm_b_threshold = 2;
sequence_param_struct.bound_iterate = 1;
sequence_param_struct.clip_threshold = 10;
sequence_param_struct.function_eval_per_iteration = 4 + ...
    sequence_param_struct.use_greedy_algorithm_b;
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
% It takes 77m to run this script. < 2Hr
% I need to fix the convergence of the algorithms.
for budget=[2000 10000];
    % Set A to be 10% of the number of iterations performed.
    sequence_param_struct.A = ...
        (budget / sequence_param_struct.function_eval_per_iteration) / 10;
for p=10:10:100 % for multiple dimensions.
    init_theta = 0.2 * ones(p, 1);
    true_loss_fn = quartic_loss_factory(p);
    target_fn = noisy_function_factory(true_loss_fn, sigma);
    true_optimal_theta = zeros(p, 1);
    for run_idx=1:50 % for multiple runs.
        for name_fn_idx=1:length(name_fn_cell) % for different algorithms
            name_fn = name_fn_cell{name_fn_idx};
            FN = name_fn_struct.(name_fn);
            common_prefix = concat_all(name_fn, p, run_idx, budget);
            disp(common_prefix);
            % Run the algorithm.
            [iteration_count, theta, time_taken, loss_sequence, mad_sequence] = FN(...
                budget, target_fn, init_theta, true_loss_fn, ...
                true_optimal_theta, sequence_param_struct);
            % collate result.
            results_struct.([common_prefix, '_time_taken']) = time_taken;
            results_struct.([common_prefix, '_final_loss']) = ...
                loss_sequence(length(loss_sequence));
            fprintf(1, 'final_loss %f\n', results_struct.([common_prefix, '_final_loss']));
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
    save([dir_prefix '/sso_project_intermediate.mat'], 'results_struct');
end
end
save([dir_prefix '/sso_project.mat'], 'results_struct');
exit;
