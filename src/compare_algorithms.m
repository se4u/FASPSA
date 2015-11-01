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
rand('seed',31415927);
randn('seed',3111113);
sigma = 0.05;
sequence_param_cell.alpha = 0.602;
sequence_param_cell.gamma = 0.101;
sequence_param_cell.a_numerator = 1e-5;
sequence_param_cell.c_numerator = 1e-1;
budget=2000;
name_fn_struct.Adaptive2SPSA = Adaptive2SPSA;
name_fn_struct.FeedbackAdaptive2SPSA = FeedbackAdaptive2SPSA;
name_fn_struct.EfficientAdaptive2SPSA = EfficientAdaptive2SPSA;
name_fn_struct.EfficientFeedbackAdaptive2SPSA = ...
    EfficientFeedbackAdaptive2SPSA;
name_fn_cell = fieldnames(name_fn_struct);
results_struct = cell()
% for multiple dimensions.
% for multiple runs.
% for different algorithms
% Run algorithm with dimensions.
% collate result.
for p=[10 20 30 40 50 60] % for multiple dimensions.
    init_theta = 0.2 * ones(p, 1);
    true_loss_fn = quartic_loss_factory(p);
    target_fn = noisy_function_factory(true_loss_fn, sigma);
    true_optimal_theta = zeros(p, 1);
    for run_idx=1:50 % for multiple runs.
        for name_fn=name_fn_cell % for different algorithms
            name_fn = name_fn{1};
            FN = getfield(name_fn_cell, name_fn);
            common_prefix = concat_all(name_fn, p, run_idx);
            result_setter = @(key, value) setfield(...
                results_struct, [common_prefix, key], value);
            % Run the algorithm.
            [iteration_count, theta, time_taken, loss_sequence, mad_sequence] ...
                = FN(budget, target_fn, init_theta, true_loss_fn, ...
                       true_optimal_theta, sequence_param_cell);
            % collate result.
            result_setter('_time_taken', ...
                          time_taken);
            result_setter('_final_loss', ...
                          loss_sequence(length(loss_sequence)));
            result_setter('_final_mad', ...
                          mad_sequence(length(mad_sequence))); ...
            result_setter('_iteration_count', ...
                          iteration_count);
        end
    end
end
save('sso_project', name_fn_cell);
