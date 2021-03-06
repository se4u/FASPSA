function [iteration_count, theta, time_taken, loss_sequence, sqdist_sequence] = ...
    FeedbackAdaptive2SPSA(budget, target_fn, init_theta, true_loss_fn, ...
                          true_optimal_theta, sequence_param_struct)
%{
Filename    : FeedbackAdaptive2SPSA.m
Description : Improved Adaptive 2SPSA with Feedback and Weighting Mechanisms.
Author      : Pushpendre Rastogi
Created     : Wed Oct 28 14:40:57 2015 (-0400)
Last-Updated: .
By: .
Update #: 0

The paper "Feedback and Weighting Mechanisms for Improving Jacobian
Estimates in the Adaptive Simultaneous Perturbation Algorithm", James Spall,
IEEE Transactions on Automatic Control, Vol 54, No 6, (2009)
DOI:10.1109/TAC.2009.2019793 introduced enhancements to the Adaptive SPSA
algorithm for improving its convergence.

The code below implements ``2SPSA with Feedback'' which requires
only four noisy measurements of the loss function at each
iteration to estimate both the gradient and Hessian of the loss function.
See `www.jhuapl.edu/SPSA` and the paper for more details.

Note that the following implementation is very basic and does not do
any per iteration gradient/hessian averaging, initialization of 2SPSA by
running 1SPSA, parameter averaging, per iteration local search. We use the
'optimal' values for A, alpha, gamma  described in
Spall, IEEE Transactions on Aerospace and Electronic Systems, 1998, 817–823
and that are also presented in
`www.jhuapl.edu/SPSA/PDF-SPSA/Matlab-Auto_gain_sel-b-w.pdf`. The rest of the
hyper-parameters are chosen through tuning.

Inputs
======
budget : ASPSA with feedback uses 4 function evaluations per iteration to
  estimate the gradients and hessians of a function when only noisy
  evaluations are available. The budget parameter specifies the total number
  of function evaluations that would be done throughout an optimization run.

target_fn :  The ``NOISY'' function that we want to ``MINIMIZE''. In
  stochastic approximation we assume that the function evaluations
  themselves are noisy.

init_theta : The initial guess of the parameter values. The sequence of
  parameter iterates starts from this point.

%-- Arguments needed only for comparative purposes. ---------------------%
% The following arguments are used only to facilitate comparative
% analysis. Ideally one could just store all the iterates generated in
% memory and compute the outputs `loss_sequence` and `sqdist_sequence` but depending on the number of
% iterations and the dimensionality of the problem it may not be possible
% to store the sequence of parameter iterates that are generated.
true_loss_fn : The true loss function. In reality we would never know this
  but for sake of comparison we can use a known loss function and then add
  noise to it.

true_optimal_theta : The true optimal theta value that minimizes the true
  loss. This parameter only makes sense for those true_loss_fn that have a
  unique global minimizer.

sequence_param_struct : A cell with special values. See `spsa_setup.m`
Outputs
=======
iteration_count : The number of iterations of Feedback Weighted 2SPSA.
  In our implementation this value is not that useful since we always run
  to the complete budget. We don't do `early stopping` or `convergence checks`.

theta : The estimate of the optimal value of theta at the end of the budget.
  We simply report the value of `theta` at the termination of the sequence
  although a more sophisticated implementation could be to report the
  average of parameters produced during iterations.

time_taken : The time taken for the optimization.

loss_sequence : A `iteration_count + 1` length sequence that contains the true
  loss value of the sequence of iterates that was produced during
  optimization. The first value if the loss of the initial value of theta.

sqdist_sequence : A `iteration count + 1` length sequence that contains the `mean
  absolute difference` between the global optimal parameters of the true
  loss function and the parameter estimate at the kth-step in the iteration.
%}

[theta_dim, max_iterations, theta, loss_sequence, sqdist_sequence, ...
 step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, init_theta, true_optimal_theta, sequence_param_struct);
cur_loss_estimate = target_fn(theta);
loss_sequence(1) = true_loss_fn(theta);
Hbar=0;
time_taken = 0;
time_preconditioning = 0;
time_blocking = 0;
time_setup = 0;
time_linsolve = 0;
time_setup_rand = 0;
time_setup_feval = 0;
if sequence_param_struct.compare_iterations
    Hbar = 10*eye(theta_dim);
Hbar_seq = zeros(length(loss_sequence) - 1, theta_dim, theta_dim);
Hbar_theta_seq = zeros(length(loss_sequence) - 1, theta_dim);
end
% settings.sum_ck_square_ck_tilda_square = 0;
% Do the actual work.
for k=0:max_iterations-1
    tic;
    [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude, setup_time_split] = ... % , sum_ccs_update] = ...
        adaptivespsa_common(k, theta, delta_fn, perturbation_size_fn, ...
                            target_fn, sequence_param_struct);
    time_setup = time_setup + toc;
    time_setup_rand = time_setup_rand + setup_time_split.time_rand;
    time_setup_feval = time_setup_feval + setup_time_split.time_feval;
    % settings.sum_ck_square_ck_tilda_square = ...
    %     sum_ccs_update.sum_ck_square_ck_tilda_square;

    % Update Hbar
    tic
    hbar_scalar_contrib = (delta_tilda_k' * Hbar * delta_k);
    Hk_hat_minus_Phi_hat_scalar = w_k * (h_k - hbar_scalar_contrib);
    Hbar_update_half = ((Hk_hat_minus_Phi_hat_scalar/2) * delta_tilda_k) * delta_k';
    Hbar = Hbar +  (Hbar_update_half + Hbar_update_half');
    if k == 0
        Hbar = adaptivespsa_common_preconditioning(Hbar, k);
    end
    % Update Theta
    time_taken = time_taken + toc;
    tic
    Hbarbar = adaptivespsa_common_preconditioning(Hbar, k);
    time_preconditioning = time_preconditioning + toc;
    tic
    tic
    solution_of_linsolve = (Hbarbar\delta_k);
    time_linsolve = time_linsolve + toc;
    proposed_update = (step_length_fn(k)*g_k_magnitude) * solution_of_linsolve;
    proposed_theta = theta - proposed_update;
    time_taken = time_taken + toc;
    tic
    [theta, cur_loss_estimate] = greedy_algorithm_b(...
        proposed_theta, target_fn, theta, cur_loss_estimate, ...
        sequence_param_struct);
    time_blocking = time_blocking + toc;
    my_fprintf(2, ...
            ['\n FA w_k %f h_k %f |g_k| %f Hk_hat_minus_Phi_hat_scalar %f ' ...
             'rcond(Hbar) %f max(max(abs(Hbar))) %f max(abs(proposed_update)) %f'], ...
            w_k, h_k, g_k_magnitude, Hk_hat_minus_Phi_hat_scalar, ...
            rcond(Hbar), max(max(abs(Hbar))), max(abs(proposed_update)));
    loss_sequence(k+2) = true_loss_fn(theta);
    sqdist_sequence(k+2) = sqdist(theta, true_optimal_theta);
    if sequence_param_struct.compare_iterations
        Hbar_seq(k+1, :, :) = Hbar;
        Hbar_theta_seq(k+1, :) = theta;
        save('../res/FeedbackAdaptive2SPSA_Hbar_seq.mat', ...
            'Hbar_seq', 'Hbar_theta_seq');
    end
end
iteration_count = k + 1;
timing.time_taken = time_taken;
timing.time_preconditioning = time_preconditioning;
timing.time_blocking = time_blocking;
timing.time_setup = time_setup;
timing.time_linsolve = time_linsolve;
timing.time_setup_rand = time_setup_rand;
timing.time_setup_feval = time_setup_feval;
time_taken = timing;