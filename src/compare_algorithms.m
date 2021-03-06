%{
Filename    : compare_algorithms.m
Description : Code to compare the Feedback|AdaptiveSPSA algorithms with their efficient counterparts.
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 17:29:17 2015 (-0500)
Last-Updated: .
By: .
Update #: 0
%}
% close all; clear; clc;
warning('off','MATLAB:maxNumCompThreads:Deprecated');
warning('off', 'comm:system:warnobsolete:obsoleteReplace');

if ~exist('disable_print')
    global disable_print;
    disable_print = 2; % disable_print = 2 disables printing of stderr.
end
if ~exist('compare_iterations')
    global compare_iterations;
    compare_iterations = 0;
end
if maxNumCompThreads() > 1
    disp('Start matlab with -singleCompThread flag.');
    % Check if display is available.
    if ~usejava('jvm') && ~feature('ShowFigureWindows')
        exit(1);
    else
        disp('I am letting you go but remember that you cant use this run for timing');
    end

end
dir_prefix = '../res';
if exist(dir_prefix, 'dir') ~= 7
    disp(['Couldnt find directory ' dir_prefix]);
    dir_prefix = '';
    disp(['Saving files to current directory' pwd()]);
end
% rand('seed',31415927);
% randn('seed',3111113);
global noise_std;
noise_std = 0; % Noise higher than 5e-4 we can't handle.
sequence_param_struct.noise_std = noise_std;
sequence_param_struct.n_runs = 3;
if ~exist('dim_p');
    global dim_p;
    dim_p = 120;
    disp(['Set p in script to ' num2str(dim_p)]);
else
    disp(['dim_p supplied as ' num2str(dim_p)]);
end
sequence_param_struct.dim_p = dim_p;
sequence_param_struct.compare_iterations = compare_iterations;
sequence_param_struct.alpha = 1; %.602;  % a = 1
% .101 or 0.499
% setting gamma = 0.101 makes the performance of the non-efficient code
% much better than the efficient code. Setting gamma = 0.49 improves the
% performance of the efficient code more. High gamma means quick decrease
% in perturbation size, somehow the efficient method wants the
% perturbations to decay quickly? That seems counter-intuitive, what
% happens actually is that the method simply keeps hitting the ceiling of
% the greedy_algorithm_a and never actually makes a step. So the
% algorithm is not really working at all if I set gamma = 0.49. I should
% set it to gamma = 0.101 if I want the method to actually make progress.
sequence_param_struct.gamma = 0.166667; % gamma = 0.101; % gamma = 0.166667
sequence_param_struct.weight_decay_rate = 0.501;
sequence_param_struct.weight_sequence_numerator = 1/10;
% a_numerator should be tuned. In section VIII of the 2009 paper
% professor spall said that he used a = 100.
% sequence_param_struct.a_numerator = 1 * noise_std;
sequence_param_struct.c_tilda_k_multiplier = 1; % 1.1; % 1
sequence_param_struct.use_greedy_algorithm_a = 1;
sequence_param_struct.greedy_algorithm_a_threshold = 1;
sequence_param_struct.use_greedy_algorithm_b = 0;
sequence_param_struct.greedy_algorithm_b_threshold = 100;
sequence_param_struct.bound_iterate = 1;
sequence_param_struct.clip_threshold = 10;
sequence_param_struct.function_eval_per_iteration = 4 + ...
    sequence_param_struct.use_greedy_algorithm_b;
% hacky_preconditioning just reverses the step_direction if its dot
% product with the gradient direction is less than zero.
sequence_param_struct.use_hacky_preconditioning = 0;
% eps=1 is essentially gradient descent.
% eps=0 is neton descent.
sequence_param_struct.hacky_preconditioning_eps = 0;
name_fn_struct = struct();
name_fn_struct.Adaptive2SPSA = @Adaptive2SPSA;
name_fn_struct.EfficientAdaptive2SPSA = @EfficientAdaptive2SPSA;
name_fn_struct.FeedbackAdaptive2SPSA = @FeedbackAdaptive2SPSA;
name_fn_struct.EfficientFeedbackAdaptive2SPSA = ...
     @EfficientFeedbackAdaptive2SPSA;
fn_fine_struct.Adaptive2SPSA = 'time_linsolve';
fn_fine_struct.FeedbackAdaptive2SPSA = 'time_linsolve';
fn_fine_struct.EfficientAdaptive2SPSA = 'time_rank_two_update';
fn_fine_struct.EfficientFeedbackAdaptive2SPSA = 'time_rank_two_update';
name_fn_cell = fieldnames(name_fn_struct);
results_struct = struct();
% for multiple budgets.
% for multiple dimensions.
% for multiple runs.
% for different algorithms
% run algorithm with dimensions.
% budget = 200e3
for budget=(1e2 * sequence_param_struct.function_eval_per_iteration)
    n_iter = (budget / sequence_param_struct.function_eval_per_iteration);
    % Set A to be 10% of the number of iterations performed.
    sequence_param_struct.A =  n_iter / 100; % n_iter / 10; n_iter / 100
for p=[dim_p]% for multiple dimensions.
    sequence_param_struct.a_numerator = 1;
    sequence_param_struct.c_numerator = 0.01;
    true_loss_fn = @quartic_loss_fast;
    target_fn = noisy_function_factory(true_loss_fn, noise_std);
    true_optimal_theta = zeros(p, 1);
    results_struct.sequence_param_struct = ...
                sequence_param_struct;

    for run_idx=1:sequence_param_struct.n_runs % for multiple runs.
        seed_for_this_run = randint(1,1,1e6);
        % Use random initializations instead of a fixed point.
        % init_theta = 0.2 * (2 * (rand(p, 1) > 0.5) - 1);
        init_theta = 0.2 * ones(p, 1);
        % Shuffle the order of the different algorithms that we try
        % for fair timing experiments.
        algorithm_sequence_to_try = randperm(length(name_fn_cell));
        algorithm_sequence_to_try = 1:length(name_fn_cell);
        for name_fn_idx=algorithm_sequence_to_try
            rand('seed', seed_for_this_run);
            randn('seed', seed_for_this_run);
            name_fn = name_fn_cell{name_fn_idx};
            FN = name_fn_struct.(name_fn);
            common_prefix = concat_all(name_fn, p, run_idx, budget);
            my_fprintf(1, '\n %s ', common_prefix);
            % Run the algorithm.
            [iteration_count, theta, time_taken, loss_sequence, sqdist_sequence] = FN(...
                budget, target_fn, init_theta, true_loss_fn, ...
                true_optimal_theta, sequence_param_struct);
            % collate result.
            results_struct.([common_prefix, '_time_taken']) = time_taken;
            results_struct.([common_prefix, '_final_loss']) = ...
                loss_sequence(length(loss_sequence));
            my_fprintf(...
                1, ...
                ['\nIteration Count %d Init SQDIST %f Final SQDIST ' ...
                 '%f Init loss %f Final loss %f time_setup %f ' ...
                 'time_preconditioning %f  time_blocking %f time_taken %f ' ...
                 fn_fine_struct.(name_fn) ' %f ' ...
                 'time_setup_rand %f  time_setup_feval %f ' ...
                 '\n'],  ...
                iteration_count, ...
                sqdist_sequence(1), ...
                sqdist_sequence(length(sqdist_sequence)), ...
                loss_sequence(1), ...
                loss_sequence(length(loss_sequence)),...
                time_taken.time_setup, ...
                time_taken.time_preconditioning, ...
                time_taken.time_blocking, ...
                time_taken.time_taken, ...
                time_taken.(fn_fine_struct.(name_fn)), ...
                time_taken.time_setup_rand, ...
                time_taken.time_setup_feval);
            results_struct.([common_prefix, '_final_sqdist']) = ...
                sqdist_sequence(length(sqdist_sequence));
            results_struct.([common_prefix, '_iteration_count']) = ...
                iteration_count;
            results_struct.([common_prefix, '_loss_sequence']) = ...
                loss_sequence;
            results_struct.([common_prefix, '_sqdist_sequence']) = ...
                sqdist_sequence;
        end
        % Write/Overwrite intermediate result files that can be observed separately.
        save([dir_prefix '/sso_project_' num2str(dim_p) '_intermediate.mat'], 'results_struct');
    end
end
end
save([dir_prefix '/sso_project_' num2str(dim_p) '.mat'], 'results_struct');
my_fprintf(1, 'Comparison Successfully Complete');
if ~feature('ShowFigureWindows')
    exit;
end
