function Theta = pwsdde_mon(data,u)
%PWSDDE_MON Formulate the monodromy matrix of periodic orbits in PWS-DDEs 
%exploiting the COCO 'coll' toolbox
% Input:
%   data: collocation data structure
%   u: vector of continuation variables
% Output: 
%   Theta: monodromy matrix of periodic orbit

% System parameters
p_idx = data.seg{end}.p_idx;
p = u(p_idx);

% Function definitions
sys = data.sys; % system definiton data structures
nt = sys.tau_no; % number of time delays
sig = data.sig; % solution signature
pi_ej = @(j) feval(sys.e,[],[],[],sig(j),7,1);      % incoming modes at events
f_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,1,0);       % vector field in mode mj
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);    % vector field Jacobian in mode mj (wrt x and colums of xd)
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,sig(j),2,i); % event condition Jacobian at ej (wrt x and colums of xd)
Jg_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,sig(j),5,i); % jump map Jacobian at ej (wrt x and colums of xd)

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
    [xd_cn,Wd_cn,id_cn,nd_cn,WT_cn,~] = pwsdde_delay_interp(data,u,2);
    [xd_t1,Wd_t1,id_t1,nd_t1,WT_t1,~] = pwsdde_delay_interp(data,u,3);
else
    % considering state dependent delays
    [xd_cn,Wd_cn,id_cn,nd_cn,WT_cn,~,Wx_cn] = ...
        pwsdde_sd_delay_interp(data,u,2);
    [xd_t1,Wd_t1,id_t1,nd_t1,WT_t1,~,Wx_t1] = ...
        pwsdde_sd_delay_interp(data,u,3);
end  

% Unpack time delays
tau = zeros(1,nt);
tau_type = zeros(1,nt);
for i = 1:nt
    if ~sys.sd_delay 
        tau(i) = feval(sys.tau,p,i,2); % time delays
        tau_type(i) = feval(sys.tau,p,i,1);   % delay types
    else
        tau(i) = feval(sys.tau,[],p,i,5); % maximum allowed time delays
        tau_type(i) = feval(sys.tau,[],p,i,1); % delay types
    end
end
ntau = ceil(max(tau)/sum(T)); % number of periods needed to cover tau_max

% Initialize the extended Jacobian matrix (no looping around)
Jx = cell(N,N);
for i=1:N
    for j=1:N
        Jx{i,j} = sparse(prod(mps{i}.x_shp)+mps{i}.Qnum+dim+1,...
            prod(mps{j}.xbp_shp)+1);
    end
end
if isempty(ntau)
    ntau = 1;   % for systems without delays
end
Jx = repmat(Jx,1,1,ntau+1); % extended form

% Go segment by segment to fill up the Monodromy matrix
for i=1:N
    % Current vector field mode
    mi = pi_ej(i);
    nx = mps{i}.x_shp(2);
    
    % Evaluate DDE RHS and its Jacobian
    xcni = reshape(mps{i}.W*xbp{i},dim,nx); % state on collocation nodes
    fcn = zeros(size(msh{i}.fka));          % vector field evaluations on tcn
    Jfcn = zeros(size(msh{i}.fdxka));       % vector field Jacobians
    for j = 1:nx
        xcn_ij_d = squeeze(xd_cn{i}(:,:,j)); % delayed state vectors
        fcn(:,j) = f_mj(mi,xcni(:,j),xcn_ij_d);
        Jfcn(:,:,j) = Jf_mj(mi,xcni(:,j),xcn_ij_d,0);
    end

    % Jacobian of f_cond wrt x
    Jfcn = msh{i}.fdxka.*Jfcn;
    Jfcn = sparse(mps{i}.fdxrows, mps{i}.fdxcols, Jfcn(:));
    Jfx = (0.5*T(i)/mps{i}.NTST)*Jfcn*mps{i}.W - mps{i}.Wp;
    
    % Jacobian of f_cond wrt T0,T
    JfT = (0.5/mps{i}.NTST)*fcn(:);

    % combined Jacobians
    Jx{i,i,end}(1:end-dim-1,:) = Jx{i,i,end}(1:end-dim-1,:) ...
        + [Jfx JfT; mps{i}.Q sparse(mps{i}.Qnum,1)];
    
    % Vector field Jacobians with respect to delayed terms
    for  j = 1:nx
        row_idx = (j-1)*dim + (1:dim); % row index
        xcn_ij_d = squeeze(xd_cn{i}(:,:,j)); % delayed state variables
        
        for k = 1:nt % span all delays

            ik = id_cn{i,k}(j); % index of containing segment
            nk = nd_cn{i,k}(j); % index of containing periodic orbit instance
            Wdk = Wd_cn{i,k}{j}; % interpolation coefficients
            WdTk = reshape(WT_cn{i,k}(j,:,:,:),dim,N,ntau+1); % derivatives wrt T
            Jfxdk = Jf_mj(mi,xcni(:,j),xcn_ij_d,k); % vector field Jacobian matrix
            
            % wrt xbp
            Jx{i,ik,end-nk}(row_idx,1:end-1) = ...
                Jx{i,ik,end-nk}(row_idx,1:end-1) ...
                + (0.5*T(i)/mps{i}.NTST)*Jfxdk*Wdk;

            % wrt xcn (for state dependent delays)
            if sys.sd_delay
                Wdxk = reshape(Wx_cn{i,k}(j,:,:),dim,dim); % derivatives of tau wrt x
                Jx{i,i,end}(row_idx,1:end-1) = Jx{i,i,end}(row_idx,1:end-1) ...
                    + (0.5*T(i)/mps{i}.NTST)*Jfxdk*Wdxk*mps{i}.W(row_idx,:);
            end
                        
            % wrt T1
            for nkk = 0:ntau
                for n = 1:N
                    Jx{i,n,end-nkk}(row_idx,end) = ...
                        Jx{i,n,end-nkk}(row_idx,end) ...
                        + (0.5*T(i)/mps{i}.NTST)*Jfxdk*WdTk(:,n,end-nkk);
                end
            end

            % Extra terms for neutral delayed terms
            if tau_type(k) == 2
                Jx{i,ik,end-nk}(row_idx,end) = Jx{i,ik,end-nk}(row_idx,end) ...
                    - (0.5*T(i)/mps{i}.NTST)*Jfxdk*Wdk*xbp{ik}/T(ik);
            end
        end
    end

end

% Evaluate zero conditions and corresponding Jacobians for each boundary
for i=1:N
    row_idx = size(Jx{i,i,end},2)-dim-1+(1:dim); % row index of event map condition
    x1_d = squeeze(xd_t1{i}(:,:,1)); % delayed state at the boundary
    if i==N  
        % move last event map condition back by one instance
        Jx{i,1,end}(row_idx,mps{1}.x0_idx) = ...
            Jx{i,1,end}(row_idx,mps{1}.x0_idx) - eye(dim);
        Jx{i,i,end-1}(row_idx,mps{i}.x1_idx) = ...
            Jx{i,i,end-1}(row_idx,mps{i}.x1_idx) + Jg_ej(i,x1{i},x1_d,0);
    else
        Jx{i,i+1,end}(row_idx,mps{i+1}.x0_idx) = ...
            Jx{i,i+1,end}(row_idx,mps{i+1}.x0_idx) - eye(dim);
        Jx{i,i,end}(row_idx,mps{i}.x1_idx) = ...
            Jx{i,i,end}(row_idx,mps{i}.x1_idx) + Jg_ej(i,x1{i},x1_d,0);
    end
    Jx{i,i,end}(end,mps{i}.x1_idx) = Jx{i,i,end}(end,mps{i}.x1_idx) ...
        + Jh_ej(i,x1{i},x1_d,0);
    
    % Jacobians with respect to delayed terms

    for k = 1:nt % span all delays

        ik = id_t1{i,k}(1); % index of containing segment
        nk = nd_t1{i,k}(1); % index of containing periodic orbit instance
        Wdk = Wd_t1{i,k}{1}; % interpolation coefficients
        WdTk = reshape(WT_t1{i,k}(1,:,:,:),dim,N,ntau+1); % derivatives wrt T
        Jgxdk = Jg_ej(i,x1{i},x1_d,k); % event map Jacobian matrix
        Jhxdk = Jh_ej(i,x1{i},x1_d,k); % event condition Jacobian matrix
        
        % wrt xbp
        if i == N
            Jx{i,ik,end-nk-1}(row_idx,1:end-1) = ...
                Jx{i,ik,end-nk-1}(row_idx,1:end-1) + Jgxdk*Wdk;
        else
            Jx{i,ik,end-nk}(row_idx,1:end-1) = ...
                Jx{i,ik,end-nk}(row_idx,1:end-1) + Jgxdk*Wdk;
        end
        Jx{i,ik,end-nk}(end,1:end-1) = Jx{i,ik,end-nk}(end,1:end-1) ...
            + Jhxdk*Wdk;

        % wrt x1 (for state dependent delays)
        if sys.sd_delay
            Wdxk = reshape(Wx_t1{i,k}(1,:,:),dim,dim); % derivatives of tau wrt x
            Jx{i,i}(row_idx,mps{i}.x1_idx) = Jx{i,i}(row_idx,mps{i}.x1_idx) ...
                + Jgxdk*Wdxk;
            Jx{i,i}(end,mps{i}.x1_idx) = Jx{i,i}(end,mps{i}.x1_idx) ...
                + Jhxdk*Wdxk;
        end
                
        % wrt T1
        if i == N
            for nkk = 0:ntau-1
                for n = 1:N
                    Jx{i,n,end-nkk-1}(row_idx,end) = ...
                        Jx{i,n,end-nkk-1}(row_idx,end) ...
                        + Jgxdk*WdTk(:,n,end-nkk);
                end
            end
        else
            for nkk = 0:ntau
                for n = 1:N
                    Jx{i,n,end-nkk}(row_idx,end) = ...
                        Jx{i,n,end-nkk}(row_idx,end) ...
                        + Jgxdk*WdTk(:,n,end-nkk);
                end
            end
        end
        for nkk = 0:ntau
            for n = 1:N
                Jx{i,n,end-nkk}(end,end) = ...
                    Jx{i,n,end-nkk}(end,end) + Jhxdk*WdTk(:,n,end-nkk);
            end
        end

        % Extra terms for neutral delayed terms
        if tau_type(k) == 2
            if i == N
                Jx{i,ik,end-nk-1}(row_idx,end) = ...
                    Jx{i,ik,end-nk-1}(row_idx,end) ...
                    - Jgxdk*Wdk*xbp{ik}/T(ik);
            else
                Jx{i,ik,end-nk}(row_idx,end) = Jx{i,ik,end-nk}(row_idx,end) ...
                    - Jgxdk*Wdk*xbp{ik}/T(ik);
            end
            Jx{i,ik,end-nk}(end,end) = Jx{i,ik,end-nk}(end,end) ...
                - Jhxdk*Wdk*xbp{ik}/T(ik);
        end

    end

end

% Concatenate Jacobian matricies 
Jf = cell(1,ntau+1);
for i = 1:ntau+1
    Jf{i} = cell2mat(Jx(:,:,i));
end

% Add extra rows and columns if ntau>1
dimJf = size(Jf{end},1);
Znm = zeros((ntau-1)*dimJf,dimJf);
Enm = eye((ntau-1)*dimJf);
Ju_p = [Znm.' Jf{end}; -Enm Znm];
Ju_m = [Jf{1:ntau-1}, Jf{ntau}; Znm Enm];

% Calculate monodromy matrix
Theta = -Ju_p\Ju_m;

end

