function [xd,Wd,id,nd,WT,WP,Wx] = pwsdde_sd_delay_interp(data,u,...
    node_type,seg_ind)
%PWSDDE_SD_DELAY_INTERP Evaluate delayed terms, corresponding segment indicies, 
%and interpolation coefficients in the multi point boundary value
%problem of PWS-DDE periodic orbits, considering state dependent delays
% - without wrapping around to enable monodromy matrix formulation
% - also accounting for neutral delays
% Input:
%   data: collocation data structure
%   u: solution vector
%   node_type: select nodes to evaluate delayed state on (default 1)
%    -> 1: mesh nodes (tbp)
%    -> 2: collocation nodes (tcn)
%    -> 3: segment bondaries (tbp(end))
%    -> 4: last subintervals of the mesh tbp(end-NCOL+1:end)
%   seg_ind: select which segment to interpolate for (default all)
% Output:
%   xd: state at t-tau(k) {N}(dim,nt,Ntm)
%   Wd: lagrange interpolation coefficients (derivatives wrt xbp)
%       {N,nt}{Ntm}(dim,dim*(NCOL+1)*NTST)
%   id: index of segment containing t-tau(k) {N,nt}(Ntm)
%   nd: index of the periodic orbit where the tq is found without looping
%       around (0 current, 1,2,... past orbits) {N,nt}(Ntm)
%   WT: derivative of the querry point wrt T without looping around
%       {N,nt}(Ntm,dim,N,ntau+1)
%   WP: derivative of the querry point wrt parameters {N,nt}(Ntm,dim,Np)
%   Wx: derivative of the querry point wrt x {N,nt}(Ntm,dim,dim)

% Initialization
N = length(data.seg);
if nargin < 4
    seg_ind = 1:N; % by default interpolate for all segments
end
dim = data.seg{seg_ind(end)}.coll_seg.maps.x_shp(1);

% unpack the used mesh for each segment
mps = cell(1,N);    % collocation maps data structures
msh = cell(1,N);    % collocation mesh data structures
uidx = cell(1,N);   % indices of xbp
xbp = cell(1,N);    % state at mesh nodes
T = zeros(1,N);     % segment lengths
T0 = zeros(1,N);    % segment starting points
tm = cell(1,N);     % collocation nodes, mesh nodes or segment boundaries
xm = cell(1,N);     % state variables on tm
for i = seg_ind
    mps{i} = data.seg{i}.coll_seg.maps;
    msh{i} = data.seg{i}.coll_seg.mesh;
    uidx{i} = data.seg{i}.uidx;
    xbp{i} = u(uidx{i}(mps{i}.xbp_idx));
    T0(i) = u(uidx{i}(mps{i}.T0_idx));
    T(i) = u(uidx{i}(mps{i}.T_idx));
    switch node_type
        case 1
            tm{i} = T0(i) + T(i)*msh{i}.tbp; % mesh nodes
            xm{i} = reshape(xbp{i},mps{i}.xbp_shp);
        case 2
            tm{i} = T0(i) + T(i)*msh{i}.tcn; % collocation nodes
            xm{i} = reshape(mps{i}.W*xbp{i},mps{i}.x_shp);
        case 3
            tm{i} = T0(i) + T(i)*msh{i}.tbp(end); % segment boundaries
            xm{i} = reshape(u(uidx{i}(mps{i}.x1_idx)),dim,1);
        case 4
            NCOL = data.seg{i}.coll_seg.int.NCOL;
            tm{i} =  T0(i) + T(i)*msh{i}.tbp(end-NCOL:end); % final subintervals
            xm_temp = reshape(xbp{i},mps{i}.xbp_shp);
            xm{i} = xm_temp(:,end-NCOL:end);
    end
end
period = sum(T);

% initialize the system delays
p_idx = data.seg{end}.p_idx;
par = u(p_idx);
Np = length(par);  % number of parameters
sys = data.sys;
nt = sys.tau_no;   % number of delays
tau_type = zeros(1,nt);
tau_max = zeros(1,nt);
for i = 1:nt
    tau_type(i) = feval(sys.tau,[],par,i,1);   % delay types
    tau_max(i) = feval(sys.tau,[],par,i,5); % maximum allowed time delay values
end
ntau = ceil(max(tau_max)/period); % ntau*period is neccessary to cover tau_max

% initialize output variables
xd = cell(1,N);    % delayed terms in t-tau
Wd = cell(N,nt);   % Lagrange interpolation coefficients        
id = cell(N,nt);   % segment indicies containing x(t-tau)
nd = cell(N,nt);   % periodic orbit instace index (0,1,...ntau)
WT = cell(N,nt);   % interpolation derivatives wrt segment lengths
WP = cell(N,nt);   % interpolation derivatives wrt parameters
Wx = cell(N,nt);   % interpolation derivatives wrt state at t

% Span all segments
for i = seg_ind
    Ntm = size(tm{i},1); % number of querry points in current segment
    xd{i} = zeros(dim,nt,Ntm);

    % Span all delays
    for k = 1:nt
        Wd{i,k} = cell(1,Ntm);
        id{i,k} = zeros(1,Ntm);
        nd{i,k} = zeros(1,Ntm);
        WT{i,k} = zeros(Ntm,dim,N,ntau+1);
        WP{i,k} = zeros(Ntm,dim,Np);
        Wx{i,k} = zeros(Ntm,dim,dim);

        % Span all query points
        for j = 1:Ntm

            xj = xm{i}(:,j); % state vector
            tau = feval(sys.tau,xj,par,k,2); % time delay
            dp_tau = feval(sys.tau,xj,par,k,3); % time delay parameter Jacobian
            dx_tau = feval(sys.tau,xj,par,k,4); % time delay state Jacobian

            i0 = i; % past segment index
            n0 = 0; % past period index
            delta_i = tm{i}(j)-T0(i0);  % substract starting point
            delta_Ti = zeros(N,ntau+1); % derivatives wrt T(i)
            delta_Ti(i0,end) = delta_i/T(i0);% account for initial segment

            % Step back segment by segment until tau is covered
            while delta_i<tau
                % Cirle around segments
                if i0>1
                    i0 = i0-1;
                else
                    i0 = N;
                    n0 = n0+1;
                end   
                % Emergency switch
                if n0>ntau
                    warning('Maximum iteration number reached while searchig for t_tau');
                    break
                end 
                % Update delta and its derivatives
                delta_i = delta_i + T(i0);
                delta_Ti(i0,end-n0) = delta_Ti(i0,end-n0) + 1;
            end

            % Lagrange interpolation
            tq = (delta_i-tau)/T(i0);
            [W, Wp, Wpp] = pwsdde_coll_W(mps{i0},tq);
                
            % Evaluate derivative of t_tau wrt Ti
            dT = 1/T(i0)*delta_Ti;
            dT(i0,end-n0) = -(delta_i-tau-...
                delta_Ti(i0,end-n0)*T(i0))/T(i0)^2;

            % Fill up output variables
            id{i,k}(j) = i0;
            nd{i,k}(j) = n0;
            if tau_type(k) == 1
                % normal point delay
                xd{i}(:,k,j) = reshape(W*xbp{i0},dim,1);
                Wd{i,k}{j} = W;
                WP{i,k}(j,:,:) = -1/T(i0)*Wp*xbp{i0}*dp_tau;
                Wx{i,k}(j,:,:) = -1/T(i0)*Wp*xbp{i0}*dx_tau;
                for n = 1:ntau+1
                    WT{i,k}(j,:,:,n) = Wp*xbp{i0}*dT(:,n).';
                end
            else
                % neutral point delay
                xd{i}(:,k,j) = reshape(1/T(i0)*Wp*xbp{i0},dim,1);
                Wd{i,k}{j} = 1/T(i0)*Wp;
                WP{i,k}(j,:,:) = -1/T(i0)^2*Wpp*xbp{i0}*dp_tau;
                Wx{i,k}(j,:,:) = -1/T(i0)^2*Wpp*xbp{i0}*dx_tau;
                for n=1:ntau+1
                    WT{i,k}(j,:,:,n) = 1/T(i0)*Wpp*xbp{i0}*dT(:,n).';
                end
            end       
        end
    end
end


end

