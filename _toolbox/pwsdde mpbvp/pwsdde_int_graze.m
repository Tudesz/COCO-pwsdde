function [d, f, J]  = pwsdde_int_graze(prob,data,u,iq)
%PWSDDE_INT_GRAZE Critical value of the grazing condition marked by the
%event function with index iq
% Input:
%   prob: COCO problem structure
%   data: collocation data structure
%   u: vector of continuation variables
%   iq: index of the checked event condition
% Output: 
%   d: collocation data structure (copy of the input "data")
%   f: evaluation of the grazing condition
%   J: corresponding Jacobian matrix

% System parameters
p_idx = data.seg{end}.p_idx;
p = u(p_idx);
Np = length(p);

% Function definitions
sys = data.sys; % system definiton data structures
nt = sys.tau_no; % number of time delays
h_e = @(x,xd) feval(sys.e,x,xd,p,iq,1,0);    % event condition 
Jh_e = @(x,xd,i) feval(sys.e,x,xd,p,iq,2,i); % event condition Jacobian 
Jhp_e = @(x,xd) feval(sys.e,x,xd,p,iq,3,0);  % event condition parameter Jacobian

% Find delayed state variables and corresponding Jacobian matrices
if ~isfield(sys,'sd_delay')
    sys.sd_delay = false; % by default use fix point delays
end
if ~sys.sd_delay
    % fix point delays
    [xd_bp,Wd_bp,id_bp,~,WT_bp,WP_bp] = pwsdde_delay_interp(data,u,1);
else
    % considering state dependent delays
    [xd_bp,Wd_bp,id_bp,~,WT_bp,WP_bp,Wx_bp] = ...
        pwsdde_sd_delay_interp(data,u,1);
end

% Unpack solution and segment meshes
N = length(data.seg);
dim = data.seg{1}.coll_seg.maps.x_shp(1);
xbp = cell(1,N);    % state at mesh nodes
mps = cell(1,N);    % collocation maps data structures
msh = cell(1,N);    % collocation mesh data structures
uidx = cell(1,N);   % indices of segment variables in u
T = zeros(1,N);     % segment lengths
for i=1:N
    mps{i} = data.seg{i}.coll_seg.maps;
    msh{i} = data.seg{i}.coll_seg.mesh;
    uidx{i} = data.seg{i}.uidx;
    xbp{i} = u(uidx{i}(mps{i}.xbp_idx));
    T(i) = u(uidx{i}(mps{i}.T_idx));
end

% Unpack time delays
tau_type = zeros(1,nt);
for i = 1:nt
    if ~sys.sd_delay 
        tau_type(i) = feval(sys.tau,p,i,1);   % delay types
    else
        tau_type(i) = feval(sys.tau,[],p,i,1);   % delay types
    end
end

% Find the index of the critival value of the grazing condition
min_gr = zeros(1,N);
min_ind = zeros(1,N);
for i=1:N
    nx = mps{i}.xbp_shp(2);
    gri = NaN(1,nx);
    xbp_i = reshape(xbp{i},mps{i}.xbp_shp);
    for j = 2:nx-1
        % evaluate the event condition for all inner points
        gri(j) = h_e(xbp_i(:,j),squeeze(xd_bp{i}(:,:,j)));
    end
    [min_gr(i),min_ind(i)] = min(abs(gri),[],'linear');
end
[~,igr] = min(min_gr,[],'linear');
jgr = min_ind(igr);

% Evaluate the event condition and the corresponding jacobian matrix
xbp_igr = reshape(xbp{igr},mps{igr}.xbp_shp);
xgr = xbp_igr(:,jgr);
xgrd = squeeze(xd_bp{igr}(:,:,jgr));
gr_cond = h_e(xgr,xgrd);
Jgr_cond = cell(1,N);
for i=1:N
    Jgr_cond{i} = sparse(1,prod(mps{i}.xbp_shp)+2+Np);
end

% Fill up the Jacobian matrix

% wrt xbp
col_idx = (jgr-1)*dim + (1:dim); % row index
Jgr_cond{igr}(col_idx) = Jgr_cond{igr}(col_idx) + Jh_e(xgr,xgrd,0);

% wrt p
Jgr_cond{igr}(mps{igr}.p_idx) = Jgr_cond{igr}(mps{igr}.p_idx) ...
    + Jhp_e(xgr,xgrd);

 for k = 1:nt % span all delays

    ik = id_bp{igr,k}(jgr); % index of containing segment
    Wdk = Wd_bp{igr,k}{jgr}; % interpolation coefficients
    WdTk = reshape(sum(WT_bp{igr,k}(jgr,:,:,:),4),dim,N); % derivatives wrt T
    WdPk = reshape(WP_bp{igr,k}(jgr,:,:),dim,Np); % derivatives wrt p
    Jhxdk = Jh_e(xgr,xgrd,k); % vector field Jacobian matrix
    
    % wrt xbp
    Jgr_cond{ik}(mps{ik}.xbp_idx) = Jgr_cond{ik}(mps{ik}.xbp_idx) ...
        + Jhxdk*Wdk;
    
    % wrt p
    Jgr_cond{ik}(mps{ik}.p_idx) = Jgr_cond{ik}(mps{ik}.p_idx) ...
        + Jhxdk*WdPk;
    
    % wrt T1
    for n = 1:N
        Jgr_cond{n}(mps{n}.T_idx) = Jgr_cond{n}(mps{n}.T_idx) ...
            + Jhxdk*WdTk(:,n);
    end

    % wrt xgr (for state dependent delays)
    if sys.sd_delay
        Wdxk = reshape(Wx_bp{igr,k}(jgr,:,:),dim,dim); % derivatives of tau wrt x
        Jgr_cond{igr}(col_idx) = Jgr_cond{igr}(col_idx) + Jhxdk*Wdxk;
    end

    % Extra terms for neutral delayed terms
    if tau_type(k) == 2
        Jgr_cond{ik}(mps{ik}.T_idx) = Jgr_cond{ik}(mps{ik}.T_idx) ...
            - Jhxdk*Wdk*xbp{ik}/T(ik);
    end

end

% Return COCO compatible outputs
d = data;
f = gr_cond;
J = cell2mat(Jgr_cond);

end
