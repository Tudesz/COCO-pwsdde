function [px,py,xlab,ylab] = nsnl_p(t,x,xd,p)
%NSNL_P Periodic orbit plot functions for dry 2 DoF impact Duffing
%oscillator
% Input:
%   t: time mesh of periodic orbit
%   x: state vector on t
%   xd: delayed state vector on t
%   p: system parameter vector
% Output:
%   px: x coordinates to be plotted
%   py: y coordinates to be plotted
%   xlab: label of x axis
%   ylab: label of y axis

px = x(1,:);
py = x(2,:)- x(1,:);
xlab = '$x_1$ (1)';
ylab = '$x_2-x_1$ (1)';
end

