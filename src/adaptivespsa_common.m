function [w_k, h_k, delta_k, delta_tilda_k, g_k_magnitude] = adaptivespsa_common(...
    k, theta, delta_fn, perturbation_size_fn, target_fn)
delta_k = delta_fn();
delta_tilda_k = delta_fn();
c_k = perturbation_size_fn(k);
c_tilda_k = c_k;
% Compute the perturbations in theta.
ck_deltak = c_k*delta_k;
theta_plus = theta + ck_deltak;
theta_minus = theta - ck_deltak;

ctk_dtk = c_tilda_k*delta_tilda_k;
theta_plus_tilda = theta_plus + ctk_dtk;
theta_minus_tilda = theta_minus + ctk_dtk;

% Create gradient and hessian
del_yk_part1 = target_fn(theta_plus_tilda) - target_fn(theta_minus_tilda);
del_yk_part2 = target_fn(theta_plus) - target_fn(theta_minus);
del_yk = del_yk_part1 - del_yk_part2;
h_k = del_yk / (2 * c_k * c_tilda_k);
w_k = 1/(k+1);

g_k_magnitude = del_yk_part2 / c_k / 2;