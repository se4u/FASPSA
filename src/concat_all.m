function s = concat_all(varargin)
%{
Filename    : concat_all.m
Description : Convert arguments to strings and concatenate with '_' in between.
Author      : Pushpendre Rastogi
Created     : Sun Nov  1 17:30:42 2015 (-0500)
Last-Updated: .
By: .
Update #: 0
%}
s = num2str(varargin{1});
for i = 2:length(varargin)
    arg = varargin{i};
    s = [s, '_', num2str(arg)];
end