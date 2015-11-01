function mean_absolute_difference = mad(a, b)
% MAD - Return the mean absolute difference of input arguments.
% In this function, a and b are vectors, that represent the 
% parameter estimate and the true parameters. and the returned
% value is a scalar.
mean_absolute_difference = mean(abs(a - b));
