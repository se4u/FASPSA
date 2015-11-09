function [theta, cur_loss_estimate] = greedy_algorithm_b(...
    proposed_theta, target_fn, theta, cur_loss_estimate, ...
    sequence_param_struct)

algo_a_block = 0;
if sequence_param_struct.use_greedy_algorithm_a
    fprintf(2, ['\n norm(theta - proposed_theta) %f ' ...
                'sequence_param_struct.greedy_algorithm_a_threshold %f'], ...
            norm(theta - proposed_theta), ...
            sequence_param_struct.greedy_algorithm_a_threshold);
    algo_a_block = max(abs(theta - proposed_theta)) > ...
        sequence_param_struct.greedy_algorithm_a_threshold;
end

algo_b_block = 0;
new_loss = cur_loss_estimate;
if sequence_param_struct.use_greedy_algorithm_b
    new_loss = target_fn(proposed_theta);
    algo_b_block = (new_loss > cur_loss_estimate + sequence_param_struct.greedy_algorithm_b_threshold);
end

if sequence_param_struct.bound_iterate
    ct = sequence_param_struct.clip_threshold;
    proposed_theta(proposed_theta > ct) = ct;
    proposed_theta(proposed_theta < -ct) = -ct;
end
if ~(algo_a_block || algo_b_block)
    theta = proposed_theta;
    cur_loss_estimate = new_loss;
end

% bound_block = 0;
% if sequence_param_struct.bound_iterate
%     ct = sequence_param_struct.clip_threshold;
%     bound_block = (sum(theta > ct) > 0) || (sum(theta < -ct) > 0);
% end
% if ~(algo_a_block || algo_b_block || bound_block)
%     theta = proposed_theta;
%     cur_loss_estimate = new_loss;
% end
