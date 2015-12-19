function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
Hbar = (Hbar + Hbar')/2;
delta_eye = max([1e-6 * exp(-k/2) 1e-10]);
Hbarbar = Hbar / max(max(Hbar)) + eye(size(Hbar));

% Hbar_sq = Hbar * Hbar;
% m = max(max(Hbar_sq));
% Hbarbar = sqrt(m) * sqrtm(Hbar_sq/m + delta_eye * eye(size(Hbar)));

% if ~(all(isreal(Hbarbar)))
%     disp(eig(Hbar));
%     disp(delta_eye);
%     disp(eig(Hbar_sq));
%     disp(eig(Hbarbar));
%     pause;
%     exit(1);
% end