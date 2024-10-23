function [prob, data] = pwsdde_coco_prob(sys,orb_data,p_names,NTST)
%PWSDDE_COCO_PROB Create a COCO compatible problem structure for perioidic
%orbits of PWS-DDEs employing the 'coll' toolbox
% Input:
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%    -> sd_delay: if true, the considered delays are state dependent
%       (default false)
%   orb_data: starting periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> t: cell array of segment time meshes
%    -> x: solution guess on t0
%    -> p: default parameter vector
%   p_names: list of parameter names in a cell array
%   NTST: starting mesh segment number (default 10)
% Output:
%   prob: COCO compatible continuation problem structure
%   data: COCO compatible problem data structure

N = length(orb_data.sig);   % number of segments needed

% initialize a 'coll' data structure for all segments
prob = coco_prob();
p0 = reshape(orb_data.p,[],1);
src_data.p0 = p0;
src_data.pnames = p_names;
u0 = cell(N,1); % starting values of the state
for i=1:N

    % create an empty 'coll' segment
    src_data.t0 = orb_data.t{i};
    src_data.x0 = orb_data.x{i};
    if nargin>3
        [data.seg{i}, sol] = pwsdde_coll_init(prob,src_data,NTST);
    else
        [data.seg{i}, sol] = pwsdde_coll_init(prob,src_data);
    end
    u0{i} = sol.u0;

    % identify the location of these variables in the global problem
    if i==1
        data.seg{i}.uidx = 1:length(u0{i});
    else
        data.seg{i}.uidx = data.seg{i-1}.uidx(end) + (1:length(u0{i}));
    end
    pidx = data.seg{i}.coll_seg.maps.p_idx;
    data.seg{i}.p_idx = data.seg{i}.uidx(pidx);

end

% formulate the governing MP-BVP as a zero problem
u0 = [cat(1,u0{:})]; % initial solution vector
data.sys = sys;
data.sig = orb_data.sig;
prob = coco_add_func(prob, 'mpbvp', @pwsdde_mpbvp, data, 'zero',...
    'u0', u0, 'F+DF');

% define system parameters as continuation variables
prob = coco_add_pars(prob,'sys_pars',data.seg{end}.p_idx,p_names);

% add gluing conditions between different instances of the parameter vector
for i=2:N
    prob = coco_add_glue(prob,sprintf('coll.p.glue%i',i-1),...
        data.seg{i-1}.p_idx,data.seg{i}.p_idx);
end

% fix segment starting points according to segment lengths
f_sum = @(prob,data,u) deal(data, u(1)-sum(u(2:end)));
for i=1:N
    Tuidx = zeros(i-1,1);
    for j=1:i-1
        Tuidx(j) = data.seg{j}.uidx(data.seg{j}.coll_seg.maps.T_idx);
    end
    T0uidx = data.seg{i}.uidx(data.seg{i}.coll_seg.maps.T0_idx);
    fid = sprintf('coll.T.sum%i',i);
    prob = coco_add_func(prob,fid,f_sum,[],'zero','uidx',[T0uidx; Tuidx]);
end

% Add complementary monitor functions
uidx = 1:data.seg{end}.uidx(end);

% track norm and amplitudes of the solution
n = data.seg{end}.coll_seg.maps.x_shp(1);
Ampl_names = strsplit(sprintf('$|A_%i|$ ',1:n));
prob = coco_add_func(prob,'x_norm',@mon_po_norm,data,'regular',...
    {'$||U_n||$',Ampl_names{1:end-1}},'uidx',uidx);

% evaluate critical Floquet multiplier
prob = coco_add_func(prob,'floq',@mon_po_mucrit,data,'regular',...
    {'$\mu_c$','$|\mu_c|$'},'uidx',uidx);

% add an event function for detecting changes in stability
prob = coco_add_event(prob,'SC','special point','$|\mu_c|$','=',1);

% export orbit data
prob = coco_add_slot(prob,'orb_data',@coco_save_data,data,'save_full');

end
