function orb_mod = orb_convert(orb,data,guess,res)
%ORB_CONVERT Convert orbit state vector to a different signature by 
% manually picking new event locations
% Input:
%   orb: periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> t: cell array of segment time meshes (column vector)
%    -> x: solution guess on t0 (length(t0) x dim)
%    -> xd: delayed evaluation of x0 at t0-tau_i (length(t0) x dim x nt)
%    -> p: parameter vector
%    -> T: segment legths
%   data: corresponding collocation data structure (contains sys)
%   guess: number of events in new solution signature or predefined Ti1 (if
%       0, a single segment is considered, if multiple values are given 
%       they are employed as a first apporximate for the event times)
%   res: interpolation resolution for the output orbit (default 100)
% Output:
%   orb_mod: data structure of the new periodic orbit
%    -> sig: sloution signature (event list)
%    -> t0: cell array of segment time meshes (column vector)
%    -> x0: solution guess on t0 (length(t0) x dim)
%    -> p0: default parameter vector
%   err: MP-BVP error at its last evaluation

if nargin < 4
    res = 100; % default output resolution
end

% Function definitions
sys = data.sys;
h_ej = @(j,x,xd) feval(sys.e,x,xd,orb.p,j,1,0);  % event condition at ej

% Double the length of the orbit
ts = [cell2mat(orb.t.'); orb.t{end}(end)+cell2mat(orb.t.')].';
us = [cell2mat(orb.x) cell2mat(orb.x)];
us_tau = zeros(size(orb.xd{1},1),size(orb.xd{1},2),...
    length(orb.xd)*size(orb.xd{1},3));
for i = 1:length(orb.xd)
    ii = (i-1)*size(orb.xd{1},3)+(1:size(orb.xd{1},3));
    us_tau(:,:,ii) = orb.xd{i};
end
us_tau = repmat(us_tau,1,1,2);

% Evaluate event surfaces in all points
ev_tol = zeros(sys.event_no,length(ts));
for j = 1:sys.event_no
    for i = 1:length(ts)
        ut = us(:,i);
        utau = squeeze(us_tau(:,:,i));
        ev_tol(j,i) = h_ej(j,ut,utau);
    end
end

% Select events based on plots
if length(guess) == 1 && guess ~=0
    figure('Name','Event selection');
    plot(ts,ev_tol); xlabel('$t$'); ylabel('$h_i$'); hold on 
    yline(0,':k'); hold off
    ev_in = ginput(guess);
    ev_t = ev_in(:,1);
    disp(ev_t.');
elseif guess == 0
    ev_t = [0 orb.t{end}(end)];
else
    ev_t = guess; % predefined segment lengths
end
N = length(ev_t)-1;  % number of new segments

% Interpolate solution onto new mesh
[ts_u,ind_u] = unique(ts);
us_u = us(:,ind_u);
orb_mod.t = cell(1,N);
orb_mod.x = cell(1,N);
orb_mod.sig = zeros(1,N);
for j = 1:N
    [~,orb_mod.sig(j)] = min(abs(ev_tol(:,find(ts>=ev_t(j+1),1)))); % determine which event occured
    orb_mod.t{j} = linspace(ev_t(j),ev_t(j+1),res).';
    orb_mod.x{j} = zeros(length(orb_mod.t{j}),size(us_u,1));
    for i = 1:size(us_u,1)
        orb_mod.x{j}(:,i) = interp1(ts_u,us_u(i,:),orb_mod.t{j}).'; % spline interploation
    end
end
orb_mod.p = orb.p;

% Close figure
ob = findobj('type','figure','name','Event selection');
close(ob);
end


