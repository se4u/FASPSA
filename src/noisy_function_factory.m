function noisy_fn  = noisy_function_factory(loss_fn, noise_std)
%{
Filename    : noisy_function_factory.m
Description : Noises up the input loss function with additive gaussian noise.
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 16:19:10 2015 (-0500)
Last-Updated: .
By: .
Update #: 0


Inputs
======
loss_fn : The function to noise up. Vector -> Scalar
noise_std : The

Outputs
=======
noise_fn : A function that takes vector input and returns noisy scalar loss.

%}
noisy_fn = @(t) loss_fn(t) ...
    + randn(1) * noise_std ...
    + randn(1, length(t)) * t * noise_std;