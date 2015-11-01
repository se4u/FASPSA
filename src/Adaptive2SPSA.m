function [iteration_count, theta, time_taken, loss_sequence, mad_sequence] = ...
    Adaptive2SPSA(budget, target_fn, init_theta, a_numerator, ...
                  c_numerator, true_loss_fn, true_optimal_theta)
%{
Filename    : Adaptive2SPSA.m
Description : Basic Version of Adaptive 2SPSA.
Author      : Pushpendre Rastogi
Created     : Wed Oct 28 14:40:57 2015 (-0400)
Last-Updated: .
By: .
Update #: 0

Adaptive Simultaneous Perturbation Stochastic Approximation (SPSA) is a
stochastic optimization method for model free optimization of loss functions.
Its main attraction is that it approximates the gradient/hessian much faster
for a model free function than the alternative of finite differentiation.

Specifically the code below implements ``2SPSA'' which requires
only four noisy measurements of the loss function at each
iteration to estimate both the gradient and Hessian of the loss function.
See `www.jhuapl.edu/SPSA` and chapter 7 of
"Introduction to Stochastic Search and Optimization" for more details of
SPSA, Adaptive SPSA and the 2SPSA variant of Adaptive SPSA.

Also note that we use the 'optimal' values for A, alpha, gamma  described in
Spall, IEEE Transactions on Aerospace and Electronic Systems, 1998, 817–823
and that are also presented in
`www.jhuapl.edu/SPSA/PDF-SPSA/Matlab-Auto_gain_sel-b-w.pdf`. The rest of the
hyper-parameters are chosen through tuning.

Also note that the following implementation is very basic and does not do
any per iteration gradient/hessian averaging, initialization of 2SPSA by
running 1SPSA, parameter averaging, per iteration local search.

Inputs
======
budget : 2SPSA uses 4 function evaluations per iteration to estimate the
  gradients and hessians. The budget parameter specifies the total number
  of function evaluations that would be done throughout an optimization run.

target_fn :  The ``NOISY'' function that we want to ``MINIMIZE''. In
  stochastic approximation we assume that the function evaluations
  themselves are noisy.

init_theta : The initial guess of the parameter values. The sequence of
  parameter iterates starts from this point.

a_numerator : The constant used in the numerator of the sequence of step
  sizes used in the iteration. This constant is chosen through tuning and
  kept fixed after tuning. The higher this constant is, the longer our step
  sizes can be.

  a_numerator may be given a `nan` value in which case it would be replaced
  by a default of 1e-5.

c_numerator : The constant used in the numerator of the sequence of
  perturbations that are used for gradient approximations. This constant is
  also chosen through tuning and then fixed. The higher this constant is,
  the larger perturbations we are allowed to make.

  c_numerator may be given a `nan` value in which case it would be replaced
  by a default of 1e-1.
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

Outputs
=======
iteration_count : The number of iterations of 2SPSA. In this implementation
  this value is not that useful since we always run to the complete budget
  and never do `early stopping` or `convergence checks`.

theta : The estimate of the optimal value of theta at the end of the budget.
  We simply report the value of `theta` at the termination of the sequence
  although a more sophisticated implementation could be to report the
  average of parameters produced during iterations.

time_taken : The time taken for the optimization. This time should be
  iteration_count * 4 * time taken for target function evaluation plus the
  time taken for the solving a linear system of equations at each step.

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
               true_optimal_theta, a_numerator, c_numerator, 1, 0.1667);

Hbar=0;

% Do the actual work.
for k=0:max_iterations
    tic;
    [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude] = adaptivespsa_common(...
        k, theta, delta_fn, perturbation_size_fn, target_fn);
    % Update Hbar
    Hbar = (1 - w_k) * Hbar + (w_k * h_k) * symmetric(delta_tilda_k*delta_k');

    % Update Theta % This step can be made faster.
    Hbarbar = adaptivespsa_common_preconditioning(Hbar, k);
    theta = theta - (step_length_fn(k)*g_k_magnitude) * (Hbarbar\delta_k);

    time_taken = time_taken + toc;

    loss_sequence(k+2) = true_loss_fn(theta);
    mad_sequence(k+2) = mad(theta, true_optimal_theta);
end