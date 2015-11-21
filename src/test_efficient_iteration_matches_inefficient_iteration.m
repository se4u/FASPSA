clc; clear;
global compare_iterations;
compare_iterations = 1
compare_algorithms;
load ../res/Adaptive2SPSA_Hbar_seq.mat; % Hbar_seq
load ../res/EfficientAdaptive2SPSA_Bbar_seq.mat; % Bbar_seq

p = size(Hbar_seq, 2);
mm = [];
for i=1:size(Hbar_seq, 1)
    Hbar = squeeze(Hbar_seq(i, :, :));
    Bbar = squeeze(Bbar_seq(i, :, :));
    n1 = norm(Hbar * Bbar - eye(p));
    n2 = norm(Bbar * Hbar - eye(p));    
    mm = [mm, max(n1, n2)];
    my_fprintf(2, '\n %d norm(Hbar * Bbar - eye(p)) %.4g ', i, n1);
    my_fprintf(2, 'norm(Bbar * Hbar - eye(p)) %.4g \n ', n2);
end
mean(mm)
median(mm)
max(mm)

clear;
disp('To actually run this first comment out the code in adaptivespsa_preconditioning');
load ../res/FeedbackAdaptive2SPSA_Hbar_seq.mat; % Hbar_seq
load ../res/EfficientFeedbackAdaptive2SPSA_Bbar_seq.mat; % Bbar_seq

p = size(Hbar_seq, 2);
mm = [];
for i=1:size(Hbar_seq, 1)
    Hbar = squeeze(Hbar_seq(i, :, :));
    Bbar = squeeze(Bbar_seq(i, :, :));
    n1 = norm(Hbar * Bbar -eye(p));
    n2 = norm(Bbar * Hbar -eye(p));    
    mm = [mm, max(n1, n2)];
    fprintf(2, '\n %d norm(Hbar * Bbar - eye(p)) %.4g ', i, n1);
    fprintf(2, 'norm(Bbar * Hbar - eye(p)) %.4g \n ', n2);
end
mean(mm)
median(mm)
max(mm)

