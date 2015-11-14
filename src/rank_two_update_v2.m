function Bbar = rank_two_update_v2(Bbar, a, b, delta_tilda_k, delta_k)
% Stage 1 update
tmp_stage = Bbar * delta_tilda_k;
tmp_deno = a + b * (delta_k' * tmp_stage);
tmp_stage_b = Bbar' * delta_k;
Bbar = Bbar + (-b/tmp_deno) * tmp_stage * tmp_stage_b';

% Stage 2 update
tmp_stage = Bbar * delta_k;
tmp_deno = a + b * (delta_tilda_k' * tmp_stage);
tmp_stage_b = Bbar' * delta_tilda_k;
Bbar = Bbar + (-b/tmp_deno) * tmp_stage * tmp_stage_b';
