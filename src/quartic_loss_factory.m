function loss_fn = quartic_loss_factory(p, varargin)
%{
Filename    : quartic_loss_factory.m
Description : A factory method for creating the quartic loss function.
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 16:11:19 2015 (-0500)
Last-Updated: .
By: .
Update #: 0


Inputs
======
p : The dimensionality of the loss function.

varargin : They translate to the following.
B (optional) : The psd matrix parameterizing the quartic loss function.

Outputs
=======
loss_fn : This loss_fn is a function that takes a p dimensional
  real valued function and returns a scalar value.
%}
if length(varargin) == 0
    B = triu(ones(p))/p;
else
    B = varargin{1};
end
if size(B, 1) ~= p
    error('Exiting in quartic_loss_factory.m. The size of B != p');
    exit(1);
end
% M = (B'*B) ; % We could do this but making the quartic loss function
% efficient is not the point of the paper.
minima = rand(p,1);
loss_fn = @(t) (t-minima)'*(B'*B)*(t-minima) ...
          + 0.1 * sum((B*(t-minima)).^3) ...
          + 0.01 * sum((B * (t - minima)).^4);