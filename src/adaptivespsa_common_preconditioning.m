function Hbar = adaptivespsa_common_preconditioning(Hbar, k)
% delta_eye = 1e-8 * exp(-k) * eye(size(Hbar));
sa = eigs((Hbar+ Hbar')/2, 1, 'SA');
if sa < 0 
    regularizer = -sa + 1e-8;
    Hbar = Hbar + regularizer * eye(size(Hbar));
end
