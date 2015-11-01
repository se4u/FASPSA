function [theta_dim, max_iterations, theta, loss_sequence, mad_sequence, ...
          time_taken, step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, function_eval_per_iteration, init_theta, ...
               true_optimal_theta, sequence_param_cell)
%{
Filename    : spsa_setup.m
Description : A 'mixin' script that contains the common initialization
              related code that all the SPSA variants need.
Author      : Pushpendre Rastogi
Created     : Wed Oct 28 15:30:02 2015 (-0400)
Last-Updated: .
By: .
Update #: 0

See the help of the individual optimization routines which call this
function to understand the inputs and outputs.

sequence_param_cell : It contains the following keys.

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

alpha : The decay rate of the step length sequence.

gamma : The decay rate of the perturbation sequence.
%}
assert(mod(budget, function_eval_per_iteration)==0);
% Set various counters, loop length and other containers.
theta_dim = length(init_theta);
max_iterations = (budget-function_eval_per_iteration)/function_eval_per_iteration;
theta = init_theta;
loss_sequence = NaN(1, 2+max_iterations);
mad_sequence = NaN(1, 2+max_iterations);
loss_sequence(1) = true_loss_fn(theta);
mad_sequence(1) = mad(theta, true_optimal_theta);
time_taken = 0;

% Set the step - length and perturbation sequence.
A = min(100, floor(budget/20));
if isnan(sequence_param_cell)
    alpha = 0.602;
    gamma = .101;
    a_numerator = 1e-5;
    c_numerator = 1e-1;
else
    alpha = sequence_param_cell.alpha;
    gamma = sequence_param_cell.gamma;
    a_numerator = sequence_param_cell.a_numerator;
    c_numerator = sequence_param_cell.c_numerator;
end
step_length_fn=@(k) a_numerator /(k+1+A)^alpha;
perturbation_size_fn=@(k) c_numerator /(k+1)^gamma;
delta_fn=@() 2*(rand(theta_dim, 1) >= 0.5)-1;
