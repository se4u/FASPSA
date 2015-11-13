function Bbar_by_a = rank_two_update(Bbar_by_a, b, delta_tilda_k, delta_k)
% Stage 1 update
tmp_stage = Bbar_by_a * delta_tilda_k;
tmp_deno = 1 + b * (delta_k' * tmp_stage);
Bbar_by_a = Bbar_by_a - tmp_stage * ((delta_k' *(b/tmp_deno)) * Bbar_by_a);
% Stage 2 update
tmp_stage = Bbar_by_a * delta_k;
tmp_deno = 1 + b * (delta_tilda_k' * tmp_stage);
Bbar_by_a = Bbar_by_a - tmp_stage * ((delta_tilda_k'*(b/tmp_deno)) * Bbar_by_a);

%% The original code from the EfficientFeedbackAdaptive2SPSA.m
% stage1_temp = Bbar * delta_tilda_k;
% stage1_deno = 1 + b * (delta_k' * stage1_temp);
% Binv = Bbar - stage1_temp * ((delta_k' * Bbar) * (b/stage1_deno));
% stage2_temp = Binv * delta_k;
% stage2_deno = 1 + b * (delta_tilda_k' * stage2_temp);
% Bbar = Binv - stage2_temp * ((delta_tilda_k' * Binv) * (b/stage2_deno));