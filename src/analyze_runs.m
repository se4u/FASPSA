%% Mean Absolute Deviation
mytstat=@(a,b) (mean(a)-mean(b))/sqrt(var(a-b)/length(a));
critical_value=@(a, thresh) tinv(thresh, length(a)-1);
% Tell me if a is really better than b
two_sample_match_pair_one_side = @(a, b, thresh) ~isinf(mytstat(a,b)) && ...
    mytstat(a,b)>critical_value(a, thresh);

%% Code to get 1 sample Error Bars
t_nm1_alpha_by_2=@(dof, interval) tinv(interval + (1-interval)/2, dof);
delta = @(len, sd, cival) t_nm1_alpha_by_2(len-1, cival/100)*(sd/sqrt(len));
get_ci = @(mu, sd, len, cival) delta(len, sd, cival);
applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :));
applyToRows = @(func, matrix) arrayfun(applyToGivenRow(func, matrix), ...
    (1:size(matrix,1))', 'UniformOutput', 1);
get_ci_per_row = @(m) applyToRows(...
    @(row) get_ci(mean(row), std(row), length(row), 95), m);
get_error_bar = @(m) [mean(m, 2) get_ci_per_row(m)];
%% Comparative Plotting Code
YLabelCell={'Milliseconds per Iteration', 'Condition Number', ...
    'Normalized Mean Absolute Deviation From Optimal Theta',...
    'Normalized Absolute Difference from Optimal Loss'};
LossTitleCell={'Sum Product Loss','Squared Quartic Loss','Rosenbrock Loss'};
for loss_idx=2:2
    figure(loss_idx);
    for splt=1:4
        % 4 subplots time/iteration, cond_num_norm, theta_diff_norm,
        % l_diff_norm
        subplot(2, 2, splt);
        hold on;
        ye=get_error_bar(...
            reshape(...
             squeeze(Adaptive2SPSA_results(loss_idx, :, :, splt)),...
             length(PCell), Realization_Count));
        errorbar(PCell, ye(:, 1), ye(:, 2), 'bo');
        ye=get_error_bar(...
            reshape(...
             squeeze(EfficientA2SPSA_results(loss_idx, :, :, splt)),...
             length(PCell), Realization_Count));
        errorbar(PCell, ye(:, 1), ye(:, 2), 'rx');
        xlabel('Iterations')
        ylabel(YLabelCell{splt});
        legend('Original Adaptive 2SPSA', 'Efficient Adaptive 2SPSA');
        if splt == 1
            aa=squeeze(Adaptive2SPSA_results(loss_idx, :, :, splt));
            bb=squeeze(EfficientA2SPSA_results(loss_idx, :, :, splt));
            % Perform two-sample matched-pair one-sided t test
            % to check whether one algo is better than the other.
            for i=1:length(PCell)
                disp(PCell(i));
                disp(two_sample_match_pair_one_side(aa(i,:), bb(i,:), 0.99));
            end
        end
    end
    suptitle(LossTitleCell{loss_idx});
end