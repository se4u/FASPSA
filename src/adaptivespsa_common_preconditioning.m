function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
Hbar = (Hbar + Hbar')/2;
delta_eye = max([1e-6 * exp(-k/2) 1e-10]);
m = max(max(Hbar));
Hbar = Hbar / m;
Hbar_sq = Hbar * Hbar;
Hbarbar = m * sqrtm(Hbar_sq + delta_eye * eye(size(Hbar)));
if ~(all(isreal(Hbarbar)))
    disp(eig(Hbar));
    disp(delta_eye);
    disp(eig(Hbar_sq));
    disp(eig(Hbarbar));
    pause;
    exit(1);
end
% while not(all(isreal(Hbarbar)))
%     delta_eye = 10 * delta_eye;
%     my_fprintf(2, '\n delta_eye %g ', delta_eye);
%     Hbarbar = Hbarbar_f(delta_eye);
% end