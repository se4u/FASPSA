function [iteration_count, theta, time_taken, loss_sequence, mad_sequence] = ...
    EfficientFeedbackAdaptive2SPSA(budget, target_fn, init_theta, true_loss_fn, ...
                                   true_optimal_theta, sequence_param_cell)
%{
Filename    : EfficientFeedbackAdaptive2SPSA.m
Description : Basic Version of Adaptive 2SPSA.
Author      : Pushpendre Rastogi
Created     : .
Last-Updated: .
By: .
Update #: 0

Efficient Feedback Adaptive 2SPSA uses the fact that a matrix inversion can
be avoided in case of Feedback Weighted Adaptive 2SPSA by using the Matrix
Inversion Lemma.

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
% memory and compute the outputs `loss_sequence` and `mad_sequence` but depending on the number of
% iterations and the dimensionality of the problem it may not be possible
% to store the sequence of parameter iterates that are generated.
true_loss_fn : The true loss function. In reality we would never know this
  but for sake of comparison we can use a known loss function and then add
  noise to it.

true_optimal_theta : The true optimal theta value that minimizes the true
  loss. This parameter only makes sense for those true_loss_fn that have a
  unique global minimizer.

sequence_param_cell : A cell with special values. See `spsa_setup.m`
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

mad_sequence : A `iteration count + 1` length sequence that contains the `mean
  absolute difference` between the global optimal parameters of the true
  loss function and the parameter estimate at the kth-step in the iteration.
%}

[theta_dim, max_iterations, theta, loss_sequence, mad_sequence, ...
 time_taken, step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, 4, init_theta, ...
               true_optimal_theta, sequence_param_cell);

Bbar=0;
% Do the actual work.
for k=0:max_iterations
    tic;
    [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude] = adaptivespsa_common(...
        k, theta, delta_fn, perturbation_size_fn, target_fn);

    % Update Bbar
    tmp_1 = (h_k - delta_tilda_k' * Hbar * delta_k)
    tmp_2 = delta_k' * Bbar * delta_tilda_k
    Bbar = Bbar - (Bbar * delta_tilda_k) * (1/tmp_1 + tmp_2) * (delta_k' * Bbar)
    % Update Hbar
    Hk_hat_minus_Phi_hat_scalar = w_k * (h_k - (delta_tilda_k' * Hbar) * delta_k);
    Hbar = Hbar + Hk_hat_minus_Phi_hat_scalar * symmetric(delta_tilda_k * delta_k');

    % Update Theta.
    theta = theta - (step_length_fn(k)*g_k_magnitude) * (Bbar * delta_k);

    time_taken = time_taken + toc;

    loss_sequence(k+2) = true_loss_fn(theta);
    mad_sequence(k+2) = mad(theta, true_optimal_theta);
end