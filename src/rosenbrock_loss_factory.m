function loss_fn = rosenbrock_loss_factory(p, varargin)
%{
Filename    : rosenbrock_loss_factory.m
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
Outputs
=======
loss_fn : This loss_fn is a function that takes a p dimensional
  real valued function and returns a scalar value.
%}
loss_fn = @(t) sum(arrayfun(@(i)(100*(t(2*i)-t(2*i-1)^2)^2+(1-t(2*i-1))^2), 1:p/2))