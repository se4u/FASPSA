function [theta, cur_loss_estimate] = greedy_algorithm_b(...
    proposed_theta, target_fn, theta, cur_loss_estimate, ...
    sequence_param_struct)

if sequence_param_struct.use_greedy_algorithm_b
    new_loss = target_fn(proposed_theta);
    if new_loss < cur_loss_estimate + sequence_param_struct.greedy_algorithm_b_threshold
        % This is a promising direction
        theta = proposed_theta;
        cur_loss_estimate = new_loss;
    end
else
    theta = proposed_theta;
end

if sequence_param_struct.bound_iterate
    ct = sequence_param_struct.clip_threshold
    theta(theta > ct) = ct;
    theta(theta < -ct) = -ct;
end