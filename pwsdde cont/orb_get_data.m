function [orb, data] = orb_get_data(run_id,lab)
%ORB_GET_DATA Extract one periodic orbit data structure from the output of 
% COCO
% Input:
%   run_id: string identifier of the continuation run
%   lab: index of periodic orbit to extract (default 1)
% Output:
%   orb: periodic orbit data structure
%    -> u: vector of continuation parameters
%    -> sig: sloution signature (event list)
%    -> t: cell array of segment time meshes (column vector)
%    -> x: solution guess on t0 (length(t0) x dim)
%    -> xd: delayed evaluation of x0 at t0-tau_i (length(t0) x dim x nt)
%    -> p: parameter vector
%    -> T: segment legths
%    -> mu: Floquet multipliers
%    -> mu_crit: critical Floquet multiplier
%    -> v_crit: eigenvector corresponding to mu_crit 
%   data: COCO compatible data structure for plotting and new continuation
%       runs

if nargin<2
    lab = 1;
end

% load orbit chart
[data, ~] = coco_read_solution('orb_data',run_id, lab);
chart = coco_read_solution(run_id,lab,'chart');
u = chart.x;

% parameter vector
pidx = data.seg{end}.p_idx;
par = u(pidx);

% Unpack solution and segment meshes
N = length(data.seg);
tbp = cell(1,N);    % mesh nodes
xbp = cell(1,N);    % state at mesh nodes
T = zeros(1,N);     % segment lengths
for i=1:N
    mps = data.seg{i}.coll_seg.maps;
    msh = data.seg{i}.coll_seg.mesh;
    uidx = data.seg{i}.uidx;
    xbp{i} = reshape(u(uidx(mps.xbp_idx)),mps.xbp_shp);
    T0 = u(uidx(mps.T0_idx));
    T(i) = u(uidx(mps.T_idx));
    tbp{i} = T0 + T(i)*msh.tbp; % mesh nodes
end

% Find delayed terms (all points corresponding to tbp-tau_i)
if ~isfield(data.sys,'sd_delay') || ~data.sys.sd_delay
    xbp_d = pwsdde_delay_interp(data,u,1);
else
    xbp_d = pwsdde_sd_delay_interp(data,u,1); 
end

% Formulate the corresponding monodromy matrix
Th = pwsdde_mon(data,u);
[V,D] = eig(full(Th));
[~,ind] = max(abs(diag(D)));

% Output data structure
orb.u = u(1:data.seg{end}.uidx(end));
orb.sig = data.sig;
orb.t = tbp;
orb.x = xbp;
orb.xd = xbp_d;
orb.p = par;
orb.T = T;
orb.mu = diag(D);
orb.mu_crit = orb.mu(ind);
orb.v_crit = V(:,ind);

end

