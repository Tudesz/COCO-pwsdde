function [d, f, J]  = pwsdde_mpbvp(prob,data,u)
%PWSDDE_MPBVP Formulate the multi-point boundary value problem of periodic
%orbits in PWS-DDEs in a discret form with its corresponding Jacobian,
%exploiting the COCO 'coll' toolbox
% Input:
%   prob: COCO problem structure
%   data: collocation data structure
%   u: vector of continuation variables
% Output: 
%   d: collocation data structure (copy of the input "data")
%   f: evaluation of the governing system of NAEs
%   J: corresponding Jacobian matrix

% System parameters
p_idx = data.seg{end}.p_idx;
p = u(p_idx);
Np = length(p);

% Function definitions
sys = data.sys; % system definiton data structures
nt = sys.tau_no; % number of time delays
sig = data.sig; % solution signature
pi_ej = @(j) feval(sys.e,[],[],[],sig(j),7,1);      % incoming modes at events
f_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,1,0);       % vector field in mode mj
h_ej = @(j,x,xd) feval(sys.e,x,xd,p,sig(j),1,0);    % event condition at ej
g_ej = @(j,x,xd) feval(sys.e,x,xd,p,sig(j),4,0);    % jump map at ej
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);    % vector field Jacobian in mode mj (wrt x and colums of xd)
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,sig(j),2,i); % event condition Jacobian at ej (wrt x and colums of xd)
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,sig(j),5,i); % jump map Jacobian at ej (wrt x and colums of xd)
Jfp_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,3,0);     % vector field parameter Jacobian in mode mj
Jhp_ej = @(j,x,xd) feval(sys.e,x,xd,p,sig(j),3,0);  % event condition parameter Jacobian at ej
Jgp_ej = @(j,x,xd) feval(sys.e,x,xd,p,sig(j),6,0);  % jump map parameter Jacobian at ej

% Unpack solution and segment meshes
N = length(data.seg);
dim = data.seg{1}.coll_seg.maps.x_shp(1);
xbp = cell(1,N);    % state at mesh nodes
x0 = cell(1,N);     % state at segment starting points
x1 = cell(1,N);     % state at segment end points
mps = cell(1,N);    % collocation maps data structures
msh = cell(1,N);    % collocation mesh data structures
uidx = cell(1,N);   % indices of segment variables in u
T = zeros(1,N);     % segment lengths
for i=1:N
    mps{i} = data.seg{i}.coll_seg.maps;
    msh{i} = data.seg{i}.coll_seg.mesh;
    uidx{i} = data.seg{i}.uidx;
    xbp{i} = u(uidx{i}(mps{i}.xbp_idx));
    x1{i} = u(uidx{i}(mps{i}.x1_idx));
    x0{i} = u(uidx{i}(mps{i}.x0_idx));
    T(i) = u(uidx{i}(mps{i}.T_idx));
end

% Find delayed state variables and corresponding Jacobian matrices
if ~isfield(sys,'sd_delay')
    sys.sd_delay = false; % by default use fix point delays
end
if ~sys.sd_delay
    % fix point delays
    [xd_cn,Wd_cn,id_cn,~,WT_cn,WP_cn] = pwsdde_delay_interp(data,u,2);
    [xd_t1,Wd_t1,id_t1,~,WT_t1,WP_t1] = pwsdde_delay_interp(data,u,3);
else
    % considering state dependent delays
    [xd_cn,Wd_cn,id_cn,~,WT_cn,WP_cn,Wx_cn] = ...
        pwsdde_sd_delay_interp(data,u,2);
    [xd_t1,Wd_t1,id_t1,~,WT_t1,WP_t1,Wx_t1] = ...
        pwsdde_sd_delay_interp(data,u,3);
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

% Evaluate vector field and continuity conditions segment by segment
f_cond = cell(N,1);
c_cond = cell(N,1);
Jfx_cond = cell(N,N);
Jcx_cond = cell(N,N);
for i=1:N
    for j=1:N
        Jfx_cond{i,j} = sparse(prod(mps{i}.x_shp),...
            prod(mps{j}.xbp_shp)+2+Np);
        Jcx_cond{i,j} = sparse(mps{i}.Qnum,...
            prod(mps{j}.xbp_shp)+2+Np);
    end
end

for i=1:N
    % Current vector field mode
    mi = pi_ej(i);
    nx = mps{i}.x_shp(2);
    
    % Evaluate DDE RHS and its Jacobian
    xcni = reshape(mps{i}.W*xbp{i},dim,nx); % state on collocation nodes
    fcn = zeros(size(msh{i}.fka));          % vector field evaluations on tcn
    Jfcn = zeros(size(msh{i}.fdxka));       % vector field Jacobians
    Jpfcn = zeros(size(msh{i}.fdpka));      % vector field parameter Jacobians
    for j = 1:nx
        xcn_ij_d = squeeze(xd_cn{i}(:,:,j)); % delayed state vectors
        fcn(:,j) = f_mj(mi,xcni(:,j),xcn_ij_d);
        Jfcn(:,:,j) = Jf_mj(mi,xcni(:,j),xcn_ij_d,0);
        Jpfcn(:,:,j) = Jfp_mj(mi,xcni(:,j),xcn_ij_d);
    end

    % Corresponding zero conditions
    fcn = msh{i}.fka.*fcn;
    f_cond{i} = (0.5*T(i)/mps{i}.NTST)*fcn(:) - mps{i}.Wp*xbp{i};

    % Continuity conditions
    c_cond{i} = mps{i}.Q*xbp{i};

    % Jacobian of f_cond wrt x
    Jfcn = msh{i}.fdxka.*Jfcn;
    Jfcn = sparse(mps{i}.fdxrows, mps{i}.fdxcols, Jfcn(:));
    Jfx = (0.5*T(i)/mps{i}.NTST)*Jfcn*mps{i}.W - mps{i}.Wp;
    
    % Jacobian of f_cond wrt p
    Jpfcn = msh{i}.fdpka.*Jpfcn;
    Jpfcn = sparse(mps{i}.fdprows, mps{i}.fdpcols, Jpfcn(:));
    Jfp = (0.5*T(i)/mps{i}.NTST)*Jpfcn;

    % Jacobian of f_cond wrt T0,T
    JfT = (0.5/mps{i}.NTST)*fcn(:);
    JfT0 = sparse(size(JfT,1),size(JfT,2));
    
    % combined Jacobians
    Jfx_cond{i,i} = Jfx_cond{i,i}+[Jfx JfT0 JfT Jfp];
    Jcx_cond{i,i} = [mps{i}.Q  sparse(mps{i}.Qnum,2+Np)];
    
    % Vector field Jacobians with respect to delayed terms

    for  j = 1:nx
        row_idx = (j-1)*dim + (1:dim); % row index
        xcn_ij_d = squeeze(xd_cn{i}(:,:,j)); % delayed state variables
        
        for k = 1:nt % span all delays

            ik = id_cn{i,k}(j); % index of containing segment
            Wdk = Wd_cn{i,k}{j}; % interpolation coefficients
            WdTk = reshape(sum(WT_cn{i,k}(j,:,:,:),4),dim,N); % derivatives wrt T
            WdPk = reshape(WP_cn{i,k}(j,:,:),dim,Np); % derivatives wrt p
            Jfxdk = Jf_mj(mi,xcni(:,j),xcn_ij_d,k); % vector field Jacobian matrix
            
            % wrt xbp
            Jfx_cond{i,ik}(row_idx,mps{ik}.xbp_idx) = ...
                Jfx_cond{i,ik}(row_idx,mps{ik}.xbp_idx) ...
                + (0.5*T(i)/mps{i}.NTST)*Jfxdk*Wdk;
            
            % wrt p
            Jfx_cond{i,i}(row_idx, mps{ik}.p_idx) = ...
                Jfx_cond{i,i}(row_idx, mps{ik}.p_idx) ...
                + (0.5*T(i)/mps{i}.NTST)*Jfxdk*WdPk;

            % wrt xcn (for state dependent delays)
            if sys.sd_delay
                Wdxk = reshape(Wx_cn{i,k}(j,:,:),dim,dim); % derivatives of tau wrt x
                Jfx_cond{i,i}(row_idx, mps{i}.xbp_idx) = ...
                    Jfx_cond{i,i}(row_idx, mps{i}.xbp_idx) ...
                    + (0.5*T(i)/mps{i}.NTST)*Jfxdk*Wdxk*mps{i}.W(row_idx,:);
            end
            
            % wrt T1
            for n = 1:N
                Jfx_cond{i,n}(row_idx,mps{n}.T_idx) = ...
                    Jfx_cond{i,n}(row_idx,mps{n}.T_idx) ...
                    + (0.5*T(i)/mps{i}.NTST)*Jfxdk*WdTk(:,n);
            end

            % Extra terms for neutral delayed terms
            if tau_type(k) == 2
                Jfx_cond{i,ik}(row_idx,mps{ik}.T_idx) = ...
                    Jfx_cond{i,ik}(row_idx,mps{ik}.T_idx) ...
                    - (0.5*T(i)/mps{i}.NTST)*Jfxdk*Wdk*xbp{ik}/T(ik);
            end

        end
    end

end

% Evaluate event conditions and event maps
g_cond = cell(1,N);
h_cond = cell(1,N);
Jgx_cond = cell(N,N);
Jhx_cond = cell(N,N);
for i=1:N
    for j=1:N
        Jgx_cond{i,j} = sparse(dim,prod(mps{j}.xbp_shp)+2+Np);
        Jhx_cond{i,j} = sparse(1,prod(mps{j}.xbp_shp)+2+Np);
    end
end

% Evaluate zero conditions and corresponding Jacobians for each boundary

for i=1:N
    x1_d = squeeze(xd_t1{i}(:,:,1)); % delayed state at the boundary
    Jgx_cond{i,i}(:,mps{i}.x1_idx) = Jgx_cond{i,i}(:,mps{i}.x1_idx) ...
        + Jg_ej(i,x1{i},x1_d,0);
    Jgx_cond{i,i}(:,mps{i}.p_idx) = Jgx_cond{i,i}(:,mps{i}.p_idx) ...
        + Jgp_ej(i,x1{i},x1_d);
    if i==N
        g_cond{i} = g_ej(i,x1{i},x1_d) - x0{1};
        Jgx_cond{i,1}(:,mps{1}.x0_idx) = Jgx_cond{i,1}(:,mps{1}.x0_idx) ...
            -eye(dim);
    else
        g_cond{i} = g_ej(i,x1{i},x1_d) - x0{i+1};
        Jgx_cond{i,i+1}(:,mps{i+1}.x0_idx) = ...
            Jgx_cond{i,i+1}(:,mps{i+1}.x0_idx) - eye(dim);
    end
    h_cond{i} = h_ej(i,x1{i},x1_d);
    Jhx_cond{i,i}(mps{i}.x1_idx) = Jhx_cond{i,i}(mps{i}.x1_idx) ...
        + Jh_ej(i,x1{i},x1_d,0);
    Jhx_cond{i,i}(mps{i}.p_idx) = Jhx_cond{i,i}(mps{i}.p_idx) ...
        + Jhp_ej(i,x1{i},x1_d);
    
    % Jacobians with respect to delayed terms

    for k = 1:nt % span all delays

        ik = id_t1{i,k}(1); % index of containing segment
        Wdk = Wd_t1{i,k}{1}; % interpolation coefficients
        WdTk = reshape(sum(WT_t1{i,k}(1,:,:,:),4),dim,N); % derivatives wrt T
        WdPk = reshape(WP_t1{i,k}(1,:,:),dim,Np); % derivatives wrt p

        Jgxdk = Jg_ej(i,x1{i},x1_d,k); % event map Jacobian matrix
        Jhxdk = Jh_ej(i,x1{i},x1_d,k); % event condition Jacobian matrix
        
        % wrt xbp
        Jgx_cond{i,ik}(:,mps{ik}.xbp_idx) = ...
            Jgx_cond{i,ik}(:,mps{ik}.xbp_idx) + Jgxdk*Wdk;
        Jhx_cond{i,ik}(mps{ik}.xbp_idx) = ...
            Jhx_cond{i,ik}(mps{ik}.xbp_idx) + Jhxdk*Wdk;
        
        % wrt p
        Jgx_cond{i,i}(:,mps{ik}.p_idx) = Jgx_cond{i,i}(:,mps{ik}.p_idx) ...
            + Jgxdk*WdPk;
        Jhx_cond{i,i}(mps{ik}.p_idx) = Jhx_cond{i,i}(mps{ik}.p_idx) ...
            + Jhxdk*WdPk;

        % wrt x1 (for state dependent delays)
        if sys.sd_delay
            Wdxk = reshape(Wx_t1{i,k}(1,:,:),dim,dim); % derivatives of tau wrt x
            Jgx_cond{i,i}(:,mps{i}.x1_idx) = Jgx_cond{i,i}(:,mps{i}.x1_idx) ...
                + Jgxdk*Wdxk;
            Jhx_cond{i,i}(mps{i}.x1_idx) = Jhx_cond{i,i}(mps{i}.x1_idx) ...
                + Jhxdk*Wdxk;
        end
        
        % wrt T1
        for n = 1:N
            Jgx_cond{i,n}(:,mps{n}.T_idx) = Jgx_cond{i,n}(:,mps{n}.T_idx) ...
                + Jgxdk*WdTk(:,n);
            Jhx_cond{i,n}(mps{n}.T_idx) = Jhx_cond{i,n}(mps{n}.T_idx) ...
                + Jhxdk*WdTk(:,n);
        end

        % Extra terms for neutral delayed terms
        if tau_type(k) == 2
            Jgx_cond{i,ik}(:,mps{ik}.T_idx) = ...
                Jgx_cond{i,ik}(:,mps{ik}.T_idx) ...
                - Jgxdk*Wdk*xbp{ik}/T(ik);
            Jhx_cond{i,ik}(mps{ik}.T_idx) = Jhx_cond{i,ik}(mps{ik}.T_idx) ...
                - Jhxdk*Wdk*xbp{ik}/T(ik);
        end
    end


end

% Return COCO compatible outputs
d = data;
f = vertcat(f_cond{:}, c_cond{:}, g_cond{:}, h_cond{:});
J = [cell2mat(Jfx_cond); cell2mat(Jcx_cond); ...
    cell2mat(Jgx_cond); cell2mat(Jhx_cond)];

end

