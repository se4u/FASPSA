function [iteration_count, theta, time_taken, loss_sequence, sqdist_sequence] = ...
    EfficientFeedbackAdaptive2SPSA(budget, target_fn, init_theta, true_loss_fn, ...
                                   true_optimal_theta, sequence_param_struct)
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

sqdist_sequence : A `iteration count + 1` length sequence that contains the
  `mean absolute difference` between the global optimal parameters of the
  true loss function and the parameter estimate at the kth-step in the
  iteration.

% Glossary for FLOPS analysis:      (Multi, Addition)
% MVM  = Matrix * Vector            (p^2 , p^2)
% VTVM = Vector' * Vector           (p,    p)
% VSD  = Vector Scalar Division     (p,    0)
% MMS  = MMA = Matrix Matrix Substraction/Addition (0, p^2)
% VVOP = Vector Vector Outer Product (p^2, 0)
% MSD  = MAtrix Scalar Division      (P^2, 0)
% SSD/SSM = Scalar Scalar Multiplic  (1,   0)
% SSA/SSA = Scalar Scalar Addition   (0,   1)
% So the total cost of
4 MVM + 2 VTVM + 2 VSD + 3 SSD + 2 SSA + 3 MMS + 1 VVOP + 1 MSD
= (5p^2 + 4p + 3, 5p^2 + 2p + 2)
%}

[theta_dim, max_iterations, theta, loss_sequence, sqdist_sequence, ...
 step_length_fn, perturbation_size_fn, delta_fn] = ...
    spsa_setup(budget, init_theta, ...
               true_optimal_theta, sequence_param_struct);
cur_loss_estimate = target_fn(theta);
loss_sequence(1) = true_loss_fn(theta);
Bbar=0;
Hbar=0;
time_taken = 0;
time_preconditioning = 0;
time_blocking = 0;
time_setup = 0;
bbar_max = 1;
% settings.sum_ck_square_ck_tilda_square = 0;
% Do the actual work.
for k=0:max_iterations-1
    tic;
    [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude] = ... %, sum_ccs_update] = ...
        adaptivespsa_common(k, theta, delta_fn, perturbation_size_fn, ...
                            target_fn, ...
                            sequence_param_struct);
    time_setup = time_setup + toc;
    % settings.sum_ck_square_ck_tilda_square = ...
    %     sum_ccs_update.sum_ck_square_ck_tilda_square;
    %% Update Bbar
    % 1 MVM + 1 VTVM + 1 SSM + 1 SSA + 1 MVM
    % The true update is
    % Hbar = Hbar + (Hk_hat_minus_Phi_hat_scalar/2) * (uv' + vu')
    %      = Hbar + b uv' + b vu'
    % Where : b = (Hk_hat_minus_Phi_hat_scalar/2)
    %         u = delta_tilda_k; v = delta_k
    % This implies we can expedite the inverse calculation by two
    % applications of the MIL.
    % The update to Bbar requires the current value of Hbar therefore
    % compute the new Bbar before updating Hbar to save memory.
    % Let B = Hbar + buv'
    % Then B^-1 = Hbar^-1 - (Hbar^-1 buv' Hbar^-1)/(1+bv' Hbar^-1 u)
    % And  Hbar = B + b vu'
    % Then Hbar^-1 = B^-1 - (B^-1 bvu' B^-1) / (1 + b u' B^-1 v)
    tic
    Hk_hat_minus_Phi_hat_scalar = w_k * (h_k - delta_tilda_k' * Hbar * delta_k);
    b = Hk_hat_minus_Phi_hat_scalar / 2;
    Hbar_update_half = (b * delta_tilda_k) * delta_k';
    Hbar = Hbar + (Hbar_update_half + Hbar_update_half');
    if k == 0
        Hbar = adaptivespsa_common_preconditioning(Hbar, k);
        Bbar = inv(Hbar);
    else
        Bbar = rank_two_update_v2(Bbar, 1/bbar_max, b, delta_tilda_k, delta_k);
    end
    % TODO: Fix the FLOPS analysis.
    %% Update Theta.
    time_taken = time_taken + toc;
    tic
    tmp_bbar_max = max(max(abs(Bbar)));
    Bbar = Bbar / tmp_bbar_max;
    bbar_max = bbar_max * tmp_bbar_max;
    cond_bbar = adaptivespsa_common_preconditioning(Bbar, k);
    time_preconditioning = time_preconditioning + toc;
    tic
    proposed_update = (step_length_fn(k) * g_k_magnitude * bbar_max) * (cond_bbar * delta_k);
    proposed_theta = theta - proposed_update;
    time_taken = time_taken + toc;
    tic
    [theta, cur_loss_estimate] = greedy_algorithm_b(...
        proposed_theta, target_fn, theta, cur_loss_estimate, ...
        sequence_param_struct);
    time_blocking = time_blocking + toc;
    my_fprintf(2, ...
            ['\n EFA w_k %f h_k %f |g_k| %f Hk_hat_minus_Phi_hat_scalar %f ' ...
             'rcond(Bbar) %f max(abs(Bbar)) %f max(abs(proposed_update)) %f'], ...
            w_k, h_k, g_k_magnitude, Hk_hat_minus_Phi_hat_scalar, ...
            rcond(Bbar), max(abs(Bbar)), max(abs(proposed_update)));
    loss_sequence(k+2) = true_loss_fn(theta);
    sqdist_sequence(k+2) = sqdist(theta, true_optimal_theta);
end
iteration_count = k + 1;
timing.time_taken = time_taken;
timing.time_preconditioning = time_preconditioning;
timing.time_blocking = time_blocking;
timing.time_setup = time_setup;
time_taken = timing;