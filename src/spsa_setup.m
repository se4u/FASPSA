function [theta_dim, max_iterations, theta, loss_sequence, sqdist_sequence, ...
          step_length_fn, perturbation_size_fn, delta_fn]  =  spsa_setup( ...
              budget, init_theta, true_optimal_theta, sequence_param_struct)
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

sequence_param_struct : It contains the following keys.

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
fepi = sequence_param_struct.function_eval_per_iteration;
assert(mod(budget, fepi)==0);
% Set various counters, loop length and other containers.
theta_dim = length(init_theta);
max_iterations = budget/fepi;
theta = init_theta;
loss_sequence = NaN(1, 1 + max_iterations);
sqdist_sequence = NaN(size(loss_sequence));
% loss_sequence(1) = true_loss_fn(theta);
sqdist_sequence(1) = sqdist(theta, true_optimal_theta);

% Set the step - length and perturbation sequence.
% if ~isfield(sequence_param_struct, 'alpha')
%     alpha = 0.602;
%     gamma = .101;
%     a_numerator = 100;
%     c_numerator = .05;
%     A = 100;
% else
alpha = sequence_param_struct.alpha;
gamma = sequence_param_struct.gamma;
a_numerator = sequence_param_struct.a_numerator;
c_numerator = sequence_param_struct.c_numerator;
A = sequence_param_struct.A;
% end
step_length_fn=@(k) a_numerator /(k+1+A)^alpha;
perturbation_size_fn=@(k) c_numerator /(k+1)^gamma;
delta_fn=@() 2*(rand(theta_dim, 1) >= 0.5)-1;
