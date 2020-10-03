function x = ceilfix(x)
% x = ceilfix(x)
% 
% Rounds values in x to the nearest integer AWAY from zero. (Opposite of
% fix())
%
% copyright Daniel Kimmel, 2017 Apr 11 

x = ceil(abs(x)).*sign(x);