function [theta_dim, max_iterations, theta, loss_sequence, mad_sequence, ...
          time_taken, step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, function_eval_per_iteration, init_theta, ...
               true_optimal_theta, a_numerator, c_numerator, alpha, gamma)
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
if isnan(alpha) alpha = 0.602; end;
if isnan(gamma) gamma = .101; end;
if isnan(a_numerator) a_numerator = 1e-5; end;
if isnan(c_numerator) c_numerator = 1e-1; end;
step_length_fn=@(k) a_numerator /(k+1+A)^alpha;
perturbation_size_fn=@(k) c_numerator /(k+1)^gamma;
delta_fn=@() 2*(rand(theta_dim, 1) >= 0.5)-1;
