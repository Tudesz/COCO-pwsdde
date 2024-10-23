function [d, f, J]  = pwsdde_ext_graze(prob,data,u,iq)
%PWSDDE_EXT_GRAZES Grazing condition of the event at the boundary marked 
% by iq in the MPBVP of PWSDDE periodic orbits
% Input:
%   prob: COCO problem structure
%   data: collocation data structure
%   u: vector of continuation variables
%   iq: index of the checked boundary
% Output: 
%   d: collocation data structure (copy of the input "data")
%   f: evaluation of the grazing condition
%   J: corresponding Jacobian matrix

% System parameters
p_idx = data.seg{end}.p_idx;
p = u(p_idx);
Np = length(p);
N = length(data.seg);

% Function definitions
sys = data.sys; % system definiton data structures
nt = sys.tau_no; % number of time delays
sig = data.sig; % solution signature
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,sig(j),1,0);    % event condition at ej
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,sig(j),2,i); % event condition Jacobian at ej (wrt x and colums of xd)
Jhp_ej = @(j,x,xd) feval(sys.e,x,xd,p,sig(j),3,0);  % event condition parameter Jacobian at ej

% Unpack solution segment and corresponding mesh data
uidx = data.seg{iq}.uidx;
mps = data.seg{iq}.coll_seg.maps;
msh = data.seg{iq}.coll_seg.mesh;
NCOL = data.seg{iq}.coll_seg.int.NCOL;
xbp = u(uidx(mps.xbp_idx));
xbp_shp = data.seg{iq}.coll_seg.maps.xbp_shp;
dim = xbp_shp(1);
% T = u(uidx(mps.T_idx));

% Unpack time delays
if ~isfield(sys,'sd_delay')
    sys.sd_delay = false; % by default use fix point delays
end
tau_type = zeros(1,nt);
for i = 1:nt
    if ~sys.sd_delay 
        tau_type(i) = feval(sys.tau,p,i,1);   % delay types
    else
        tau_type(i) = feval(sys.tau,[],p,i,1);   % delay types
    end
end

% Find delayed state variables and corresponding Jacobian matrices
if ~sys.sd_delay
    % fix point delays
    [xd_bi,Wd_bi,id_bi,~,WT_bi,WP_bi] = pwsdde_delay_interp(data,u,4,iq);
else
    % considering state dependent delays
    [xd_bi,Wd_bi,id_bi,~,WT_bi,WP_bi,Wx_bi] = ...
        pwsdde_sd_delay_interp(data,u,4,iq);
end

% Evaluate event condition for all points in current segment
hbp = zeros(NCOL+1,1); % event conditions
xbp_vec = reshape(xbp,xbp_shp);
for i = 1:NCOL+1
    xbpi = xbp_vec(:,end-NCOL-1+i);     % state vector
    xbpi_d = squeeze(xd_bi{iq}(:,:,i));  % delayed state vectors
    hbp(i) = h_ej(iq,xbpi,xbpi_d);
end

% Formulate grazing condtion
[~, Wp, ~] = pwsdde_coll_W(mps,msh.tbp(end)); % derivative matrix for the segment boundary
col_ind = prod(xbp_shp)-dim*(NCOL+1) + 1 + dim*(0:NCOL);
Wp1 = Wp(1,col_ind); % select only relevant variables (COULD BE IMPROVED WITH A PROBLEM SPECIFIC ALGORITHM)
gr_cond = Wp1*hbp; % approximate derivative of h(tbp(end))

% Corresponding Jacobian matricies
Jgr_cond = cell(1,N);
for i=1:N
    mps_temp = data.seg{i}.coll_seg.maps;
    Jgr_cond{i} = sparse(1,prod(mps_temp.xbp_shp)+2+Np);
end

for i = 1:NCOL+1
    xbpi = xbp_vec(:,end-NCOL-1+i);     % state vector
    xbpi_d = squeeze(xd_bi{iq}(:,:,i));  % delayed state vectors

    % wrt xbp
    col_ind = prod(xbp_shp)-dim*(NCOL+1) + dim*(i-1) + (1:dim);
    Jgr_cond{iq}(col_ind) =  Jgr_cond{iq}(col_ind) ...
        + Wp1(i)*Jh_ej(iq,xbpi,xbpi_d,0);

    % wrt p
    p_ind = prod(xbp_shp) + 2 + (1:Np);
    Jgr_cond{iq}(p_ind) =  Jgr_cond{iq}(p_ind) ...
        + Wp1(i)*Jhp_ej(iq,xbpi,xbpi_d);

   for k = 1:nt % span all delays

            ik = id_bi{iq,k}(i); % index of containing segment
            Wdk = Wd_bi{iq,k}{i}; % interpolation coefficients
            WdTk = reshape(sum(WT_bi{iq,k}(i,:,:,:),4),dim,N); % derivatives wrt T
            WdPk = reshape(WP_bi{iq,k}(i,:,:),dim,Np); % derivatives wrt p
            Jhxdk = Jh_ej(iq,xbpi,xbpi_d,k); % event condition Jacobian matrix
            mps_k = data.seg{ik}.coll_seg.maps;

            % wrt xbp
            Jgr_cond{ik}(mps_k.xbp_idx) = Jgr_cond{ik}(mps_k.xbp_idx) ...
                + Wp1(i)*Jhxdk*Wdk;

            % wrt p
            Jgr_cond{ik}(mps_k.p_idx) = Jgr_cond{ik}(mps_k.p_idx) ...
                + Wp1(i)*Jhxdk*WdPk;
            
            % wrt xbp (for state dependent delays)
            if sys.sd_delay
                Wdxk = reshape(Wx_bi{iq,k}(i,:,:),dim,dim); % derivatives of tau wrt x
                Jgr_cond{iq}(col_ind) =  Jgr_cond{iq}(col_ind) ...
                    + Wp1(i)*Jhxdk*Wdxk;
             end

            % wrt T1
            for n = 1:N
                mps_n = data.seg{n}.coll_seg.maps;
                Jgr_cond{n}(mps_n.T_idx) = Jgr_cond{n}(mps_n.T_idx) ...
                    + Wp1(i)*Jhxdk*WdTk(:,n);
            end

            % Extra terms for neutral delayed terms
            if tau_type(k) == 2
                xbpk = u(data.seg{ik}.uidx(mps_k.xbp_idx));
                Tk = u(data.seg{ik}.uidx(mps_k.T_idx));
                Jgr_cond{ik}(mps_k.T_idx) = Jgr_cond{ik}(mps_k.T_idx)...
                    - Wp1(i)*Jhxdk*Wdk*xbpk/Tk;
            end

   end

end

% Return COCO compatible outputs
d = data;
f = gr_cond;
J = cell2mat(Jgr_cond);

end




