function [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude, varargout] = adaptivespsa_common(...
    k, theta, delta_fn, perturbation_size_fn, target_fn, sequence_param_struct, ...
    varargin)
%{
Filename    : adaptivespsa_common.m
Description : Common code used by all the adaptivespsa variants.
Author      : Pushpendre Rastogi
Created     : Wed Nov  4 20:05:32 2015 (-0500)
Last-Updated: .
By: .
Update #: 0

Adaptive SPSA and its variants

Inputs
======
k : the iteration index (it starts at 0 in the code, but 1 in the paper).

theta : the current estimate of the values of the parameters.

delta_fn : The function used to produce the binary vector of perturbations
  in the parameters.

perturbation_size_fn : The function that produces the value of c_k at iteration k.

target_fn : The noisy target function.

sequence_param_struct : It contains the following fields:
  c_tilda_k_multiplier : the paper "SPSA by simultaneous perturbation method"
    says that setting c_tilda > c is better. That leaves open the possibility
    that c_tilda may be set higher to something like c_tilda =  1.1 * c
    or 2 * c ?
  weight_decay_rate : A way to set the decay rate of the weights.
varargin
========
varargin{1}.sum_ck_square_ck_tilda_square
sum_ck_square_ck_tilda_square : The optimal gain sequence in the case of
  feedback weighted adaptive SPSA requires the computation of a sequence of
  partial sums of \sum_k (c_k c_tilda_k)^2. We pass in that sequence through
  the 1st argument of varargin.

Outputs
=======
w_k : The weighting parameter for updating Hbar_k

h_k : Short form for del_yk / (2 * c_k * c_tilda_k).

delta_k : The perturbation vector. A binary vector of -1, 1

delta_tilda_k : The perturbation vector. A binary vector of -1, 1

g_k_magnitude : The scalar multiplier to delta_k^-1 which is a biased
  estimate of G_k(theta_k)

varargout : return the updated sum_ck_square_ck_tilda_square.
%}
delta_k = delta_fn();
delta_tilda_k = delta_fn();
c_k = perturbation_size_fn(k);
c_tilda_k = c_k * sequence_param_struct.c_tilda_k_multiplier;

% Compute the perturbations in theta.
ck_deltak = c_k * delta_k;
theta_plus = theta + ck_deltak;
theta_minus = theta - ck_deltak;

ctk_dtk = c_tilda_k * delta_tilda_k;
theta_plus_tilda = theta_plus + ctk_dtk;
theta_minus_tilda = theta_minus + ctk_dtk;

% Create gradient and hessian
del_yk_part1 = target_fn(theta_plus_tilda) - target_fn(theta_minus_tilda);
% del_yk_part2 will be reused later to compute the gradient estimate.
del_yk_part2 = target_fn(theta_plus) - target_fn(theta_minus);
del_yk = del_yk_part1 - del_yk_part2;
% Note that if del_yk does not become small as quickly as c_k^2 then
% h_k would rise fast. If h_k rises faster than w_k can compensate
% then the matrix would become more and more ill-conditioned.
ck_ctk = c_k * c_tilda_k;
fprintf(2, '\n ck_ctk %f del_yk %f ', ck_ctk, del_yk);
h_k = del_yk / (2 * ck_ctk);

if length(varargin) == 0
    % The w_k sequence 1/(k+1) should be used for the basic Adaptive2SPSA.
    w_k = 1/(k+1);
else
    % The optimal weighting in case of feedback weighted adaptive SPSA is
    % given in professor Spall's 2009 paper in the equation (4.2)
    ck_square_ck_tilda_square = (c_k * c_tilda_k)^2;
    sum_ck_square_ck_tilda_square = varargin{1}.sum_ck_square_ck_tilda_square + ...
        ck_square_ck_tilda_square;
    varargout{1}.sum_ck_square_ck_tilda_square = sum_ck_square_ck_tilda_square;
    w_k = ck_square_ck_tilda_square / sum_ck_square_ck_tilda_square;
end
g_k_magnitude = del_yk_part2 / c_k / 2;