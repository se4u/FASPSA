%% SSO Project
close all; clear; clc;
%% Mean Absolute Deviation
MAD = @(t1, t2) mean(abs(t1-t2));
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

%% Code to Get Comparison Data points
PCell=[10 20 30 40 50 60];
budget = 8000;
Realization_Count=20;
Adaptive2SPSA_results = nan(3, length(PCell), ...
    Realization_Count, 4);
EfficientA2SPSA_results = nan(size(Adaptive2SPSA_results));
for p_idx=1:length(PCell)
    p=PCell(p_idx);
    B = triu(ones(p))/p;
    ThetaOptCell = {exp(sum(log(1:p))/p)*(1./(1:p)'), ... %SumProd 
        zeros(p, 1), ... %Quartic
        ones(p, 1)}; %Rosenbrock
    LCell = {@(t) sum((1:p).*(t')) + prod(t.^(-1)), ... %SumProd
        @(t) t'*(B'*B)*t + 0.1 * sum((B*t).^3) + 0.01 * sum((B * t).^4), ... 
        @(t) sum(arrayfun(@(i)(100*(t(2*i)-t(2*i-1)^2)^2+(1-t(2*i-1))^2), 1:p/2))};
    for loss_idx=2:2%length(LCell)
        L=LCell{loss_idx};
        OptTheta=ThetaOptCell{loss_idx};
        OptL = L(OptTheta);
        InitTheta = 1.1*OptTheta+ones(p,1);
        InitL = L(InitTheta);
        for noise_sigma=0.01
            Y = @(t) L(t) + randn(1)*noise_sigma;
            rand('seed',31415927);
            randn('seed',3111113);    
            for rlztn = 1:Realization_Count
                % Efficient 2SPSA using matrix inversion lemma
                [iterations, theta_hat, condition_number, time_taken, ls, md] = ...
                    EfficientAdaptive2SPSA(budget, Y, InitTheta, p,...
                                              1, 4, L, OptTheta, OptL);
                EfficientA2SPSA_results(loss_idx, p_idx, rlztn, :) = ...
                    [time_taken/iterations*1000, condition_number, ...
                     md(length(md)), ls(length(ls))];
            end
            rand('seed',31415927);
            randn('seed',3111113);    
            for rlztn = 1:Realization_Count
                % Original 2SPSA
                [iterations, theta_hat, condition_number, time_taken, ls, md] = ...
                    Adaptive2SPSA(budget, Y, InitTheta, p, ...
                                     1, 4, L, OptTheta, OptL);
                Adaptive2SPSA_results(loss_idx, p_idx, rlztn, :) = ...
                    [time_taken/iterations*1000, condition_number, ...
                     md(length(md)), ls(length(ls))];
            end
        end
    end
end
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