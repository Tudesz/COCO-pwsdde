function [px,py,xlab,ylab] = nlne_robot_p(t,x,xd,p)
%NLNE_ROBOT_P Periodic orbit plot functions for nonlinear robotic arm 
% Input:
%   t: time mesh of periodic orbit
%   x: state vector on t
%   xd: delayed state vector on t
%   p: system parameter vector
% Output:
%   px: x coordinates to be plotted
%   py: y coordinates to be plotted

px = x(1,:);
py = x(2,:);
xlab = '$x_1$ (1)';
ylab = '$\dot{x}_1$ (1)';

end