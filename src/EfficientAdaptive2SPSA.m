function [iteration_count, theta, time_taken, loss_sequence, sqdist_sequence] = ...
    EfficientAdaptive2SPSA(budget, target_fn, init_theta, true_loss_fn, ...
                           true_optimal_theta, sequence_param_struct)
%{
Filename    : EfficientAdaptive2SPSA.m
Description : Basic Version of Adaptive 2SPSA.
Author      : Pushpendre Rastogi
Created     : .
Last-Updated: .
By: .
Update #: 0

Efficient Adaptive 2SPSA uses the fact that a matrix inversion can be avoided in
case of Adaptive 2SPSA by using the Matrix Inversion Lemma.

Inputs
======
budget : The budget parameter specifies the total number
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

Outputs
=======
iteration_count : The number of iterations of 2SPSA. In this implementation
  this value is not that useful since we always run to the complete budget
  and never do `early stopping`, `Blocking`, or `convergence checks`.

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

sequence_param_struct : A cell with special values. See `spsa_setup.m`
%}

[theta_dim, max_iterations, theta, loss_sequence, sqdist_sequence, ...
 step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, init_theta, ...
               true_optimal_theta, sequence_param_struct);
cur_loss_estimate = target_fn(theta);
loss_sequence(1) = true_loss_fn(theta);
Bbar = 0;
time_taken = 0;
time_preconditioning = 0;
time_blocking = 0;
time_setup = 0;
time_rank_two_update = 0;
bbar_max = 1;
% Do the actual work.
for k=0:max_iterations-1
    tic;
    [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude] = adaptivespsa_common(...
        k, theta, delta_fn, perturbation_size_fn, target_fn, ...
        sequence_param_struct);
    time_setup = time_setup + toc;
    %% Update Bbar
    % The primal update is
    % Hbar = (1 - w_k) * Hbar + (w_k * h_k) * symmetric(delta_tilda_k * delta_k')
    % This is the same as:
    % Hbar = a Hbar + b uv' + b vu';
    % Where: a = 1 - w_k; b = (w_k * h_k)/2
    %        u = delta_tilda_k; v = delta_k
    % And we would implement the inversion as a series of two rank 1 updates.
    % Let  B = a Hbar + b uv'
    % B^-1  = Hbar^-1 / a  - ((Hbar^-1/a) buv' (Hbar^-1/a))/(1 + bv'(Hbar^-1/a)u)
    % Then Hbar = B + b vu'
    % Then Hbar^-1 = B^-1 - (B^-1 bvu' B^-1)/(1 + b u' B^-1 v)
    %% Prep
    tic
    a = 1 - w_k;
    b = w_k * h_k / 2;
    if k == 0
        Bbar = inv(adaptivespsa_common_preconditioning(...
            b * ((delta_tilda_k * delta_k') + (delta_k * delta_tilda_k')), k));
    else
        tic
        Bbar = rank_two_update_v2_fast(Bbar, a/bbar_max, b, delta_tilda_k, delta_k);
        time_rank_two_update = time_rank_two_update + toc;
    end
    % Update Theta.
    % It is critical to use this modified newton step of converting
    % the negative eigen values of Bbar to positive.
    % Right now we are doing an expensive operation but this can be sped
    % up considerably.
    tmp_bbar_max = max(max(abs(Bbar)));
    Bbar = Bbar / tmp_bbar_max;
    bbar_max = bbar_max * tmp_bbar_max;
    time_taken = time_taken + toc;
    tic
    cond_bbar = adaptivespsa_common_preconditioning(Bbar, k);
    time_preconditioning = time_preconditioning + toc;
    tic
    step_size = (step_length_fn(k)*g_k_magnitude*bbar_max);
    %% The No Preconditioning Solution.
    % step_direction = Bbar * delta_k;
    % if step_direction'*delta_k < 0
    %     step_size = -step_size;
    % end
    step_direction = cond_bbar * delta_k;
    proposed_theta = theta - step_size * step_direction;
    time_taken = time_taken + toc;
    tic
    [theta, cur_loss_estimate] = greedy_algorithm_b(...
        proposed_theta, target_fn, theta, cur_loss_estimate, ...
        sequence_param_struct);
    time_blocking = time_blocking + toc;
    loss_sequence(k+2) = true_loss_fn(theta);
    sqdist_sequence(k+2) = sqdist(theta, true_optimal_theta);
end
iteration_count = k + 1;
timing.time_taken = time_taken;
timing.time_preconditioning = time_preconditioning;
timing.time_blocking = time_blocking;
timing.time_setup = time_setup;
timing.time_rank_two_update = time_rank_two_update;
time_taken = timing;