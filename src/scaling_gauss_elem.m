clear; close all; clc
pv = [  5 6] * 1e4; % 1 2 4
tv = nan(length(pv), 1);
runs=1;
tv2 = nan(length(pv), 1);
for pi = 1:length(pv)
	   p = pv(pi);
v = randn(p, 1);
m = randn(p, p);
    % Calculate time for mldivide.
    % tic; for i = 1:runs mldivide(m, v); end; tv(pi) = toc/runs ;
    % Calculate time for linsolve.
    tic; for i = 1:runs linsolve(m, v); end; tv2(pi) = toc/runs;
end

% fprintf(1, '%f ', tv); fprintf(1, '\n');
fprintf(1, '%f ', tv2); fprintf(1, '\n');
%plot(pv, tv, 'b');
%hold on;
%plot(pv, tv2, 'r');
%ylabel('nano-seconds');
%xlabel('dimension');
