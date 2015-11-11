function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
delta_eye = 1e-8 * exp(-k) * eye(size(Hbar));
Hbar_sq = Hbar * Hbar;
Hbarbar = i;
Hbarbar_f  = @(regularizer) sqrtm(Hbar_sq + regularizer);
while not(all(isreal(Hbarbar)))
    Hbarbar = Hbarbar_f(delta_eye);
    delta_eye = 10 * delta_eye;
end
% % Hbarbar = Hbar + delta_eye; % This conditioning is incredibly bad!
% fprintf(2, '\n rcond(Hbar) %f rcond(Hbarbar) %f', ...
%         rcond(Hbar), rcond(Hbarbar));