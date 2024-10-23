function plot_orb(orb,fs)
%PLOT_ORB Plot the state or defined functions on the bvp solution
% Input:
%   orb: periodic orbit data structure (based on the simulation guess)
%    -> sig: sloution signature (event list)
%    -> t: cell array of segment time meshes (column vector)
%    -> x: solution guess on t0 (length(t0) x dim)
%    -> xd: delayed evaluation of x0 at t0-tau_i (length(t0) x dim x nt)
%    -> p: parameter vector
%    -> T: segment legths
%   fs: functions for plotting, default: fs1(t,x)=t wrt fs2(t,x)=x

% Unpack orbit data (working for the output of sim_ns_dde.m as well)
t = vertcat(orb.t{:}).';
nt = length(orb.t{1});
if size(orb.x{1},1) == nt
    x = vertcat(orb.x{:}).';
else
    x = horzcat(orb.x{:});
end
if isfield(orb,'xd')
    xd = cat(3,orb.xd{:});
else
    xd = zeros([size(x) 0]);
end
if isfield(orb,'p0')
    p = orb.p0;
else
    p = orb.p;
end
if isfield(orb,'T')
   T = orb.T;
else
   T = diff([0 t(nt:nt:end)]);
end

if nargin<2
    plot(t,x,'.-');
    for i=1:length(T)+1
        xline(sum(T(1:i-1)),':k');
        xlabel('$t$'); ylabel('$x$');
    end
else
    [px,py,xlab,ylab] = feval(fs,t,x,xd,p);
    plot(px,py,'k.-'); hold on;
    xlabel(xlab); ylabel(ylab);
end

hold off;
box on;

end

