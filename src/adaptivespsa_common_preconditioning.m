function Hbarbar = adaptivespsa_common_preconditioning(Hbar, k)
% delta_eye = 1e-8 * exp(-k) * eye(size(Hbar));
delta_eye = 1e-10;
Hbar_sq = Hbar * Hbar;
Hbarbar_f  = @(regularizer) sqrtm(Hbar_sq + regularizer * eye(size(Hbar)));
Hbarbar = Hbarbar_f(delta_eye);
while not(all(isreal(Hbarbar)))
    delta_eye = 10 * delta_eye;
    fprintf(2, '\n delta_eye %g ', delta_eye);
    Hbarbar = Hbarbar_f(delta_eye);
end
% % Hbarbar = Hbar + delta_eye; % This conditioning is incredibly bad!
% fprintf(2, '\n rcond(Hbar) %f rcond(Hbarbar) %f', ...
%         rcond(Hbar), rcond(Hbarbar));