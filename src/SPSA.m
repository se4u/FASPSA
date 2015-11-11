function [iteration_count, theta, time_taken, loss_sequence, sqdist_sequence] = ...
    SPSA(budget, target_fn, init_theta, true_loss_fn, ...
         true_optimal_theta, sequence_param_cell)
%{
Filename : SPSA.m
Description : First Order SPSA for model free optimization.
Author : Pushpendre Rastogi
Created : Wed Oct 28 12:30:52 2015 (-0400)
Last-Updated: .
By: .
Update #: 0

Simultaneous Perturbation Stochastic Approximation (SPSA) is a stochastic
optimization method for model free optimization of loss functions. Its main
attraction is that it approximates the gradient/hessian much faster
for a model free function than the alternative of finite differentiation.
See `www.jhuapl.edu/SPSA` for more details.

Specifically the code below implements the algorithm presented in
`www.jhuapl.edu/SPSA/PDF-SPSA/Matlab-SPSA_Alg.pdf`
And it uses the 'optimal' values for A, alpha, gamma that were described in
Spall, IEEE Transactions on Aerospace and Electronic Systems, 1998, 817–823
and that are also presented in
`www.jhuapl.edu/SPSA/PDF-SPSA/Matlab-Auto_gain_sel-b-w.pdf`. The rest of the
values are chosen through tuning.

Inputs
======
budget : Order 1 SPSA uses function evaluations to approximate the gradient
  and then uses SGD style parameter updates. The budget refers to the number
  of function evaluations done before we stop iterating and return with the
  final parameters.

target_fn : The ``NOISY'' function that we want to ``MINIMIZE''. In
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

sequence_param_cell : A cell with special values. See `spsa_setup.m`
Outputs
=======
iteration_count : The number of iterations of 1SPSA. In this implementation
  this value is not that useful since we always run to the complete budget
  and never do `early stopping` or `convergence checks`.

theta : The estimate of the optimal value of theta at the end of the budget.
  We simply report the value of `theta` at the termination of the sequence
  although a more sophisticated implementation could be to report the
  average of parameters produced during iterations.

time_taken : The time taken for the optimization. This time should just be
  iteration_count * 2 * time taken for target function evaluation.

loss_sequence : A `iteration_count + 1` length sequence that contains the true
  loss value of the sequence of iterates that was produced during
  optimization. The first value if the loss of the initial value of theta.

sqdist_sequence : A `iteration count + 1` length sequence that contains the `mean
  absolute difference` between the global optimal parameters of the true
  loss function and the parameter estimate at the kth-step in the iteration.
%}
[theta_dim, max_iterations, theta, loss_sequence, sqdist_sequence, ...
 time_taken, step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, 2, init_theta, ...
               true_optimal_theta, sequence_param_cell);

% Do the actual work.
for k=0:max_iterations-1
    tic;
    c_k = perturbation_size_fn(k);
    delta_k = delta_fn();
    tmp = c_k*delta_k;
    t_plus = theta+tmp;
    t_minus = theta-tmp;
    theta = theta - (step_length(k)/(2*c_k))*(target_fn(t_plus)-target_fn(t_minus))*delta_k;
    time_taken = time_taken + toc;
    loss_sequence(k+2)=true_loss_fn(theta);
    sqdist_sequence(k+2)=sqdist(theta, true_optimal_theta);
end
iteration_count = k + 1;
