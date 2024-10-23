function [d, f, J]  = pwsdde_slide(prob,data,u,iq)
%PWSDDE_SLIDE Sliding condition of the event at the boundary marked 
% by iq in the MPBVP of PWSDDE periodic orbits
%--------------------------------------------------------------------------
% IMPORTANT: in the current implementation, these assumptions must hold:
%   -> dh/dx is independet of x and xd
%   -> g(x) = x at the sliding event
%--------------------------------------------------------------------------
% Input:
%   prob: COCO problem structure
%   data: collocation data structure
%   u: vector of continuation variables
%   iq: index of the checked boundary
% Output: 
%   d: collocation data structure (copy of the input "data")
%   f: evaluation of the sliding condition
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
pi_ej = @(j,k) feval(sys.e,[],[],[],sig(j),7,k);    % mode change at events
f_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,1,0);       % vector field in mode mj
Jf_mj = @(mj,x,xd,i) feval(sys.f,x,xd,p,mj,2,i);    % vector field Jacobian in mode mj (wrt x and colums of xd)
Jh_ej = @(j,x,xd,i) feval(sys.e,x,xd,p,sig(j),2,i); % event condition Jacobian at ej (wrt x and colums of xd)
Jfp_mj = @(mj,x,xd) feval(sys.f,x,xd,p,mj,3,0);     % vector field parameter Jacobian in mode mj

% Unpack solution segment and corresponding mesh data
uidx = data.seg{iq}.uidx;
mps = data.seg{iq}.coll_seg.maps;
xbp = u(uidx(mps.xbp_idx));
x1 = u(uidx(mps.x1_idx));
xbp_shp = data.seg{iq}.coll_seg.maps.xbp_shp;
dim = xbp_shp(1);

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
    [xd_ti,Wd_ti,id_ti,~,WT_ti,WP_ti] = pwsdde_delay_interp(data,u,3,iq);
else
    % fix point delays
    [xd_ti,Wd_ti,id_ti,~,WT_ti,WP_ti,Wx_ti] = ...
        pwsdde_sd_delay_interp(data,u,3,iq);
end

% Evaluate sliding condition
if pi_ej(iq,1) == pi_ej(iq,2)
    warn('Sliding check without mode change, h_sl will be infinite!');
end
dh = Jh_ej(iq,x1,xd_ti{iq},0);          % event condition Jacobain
f_in = f_mj(pi_ej(iq,1),x1,xd_ti{iq});  % incoming vector field mode
f_out = f_mj(pi_ej(iq,2),x1,xd_ti{iq}); % outgoing vector field mode
h_sl = (dh*f_in + dh*f_out)/(dh*f_in - dh*f_out);

% Corresponding Jacobian matricies
Jsl_cond = cell(1,N);
for i = 1:N
    mps_temp = data.seg{i}.coll_seg.maps;
    Jsl_cond{i} = sparse(1,prod(mps_temp.xbp_shp)+2+Np);
end

% wrt xbp
Jf_in = Jf_mj(pi_ej(iq,1),x1,xd_ti{iq},0);
Jf_out = Jf_mj(pi_ej(iq,2),x1,xd_ti{iq},0);
dh_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jf_in+Jf_out))...
    -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*(dh*(Jf_in-Jf_out));
Jsl_cond{iq}(mps.x1_idx) =  Jsl_cond{iq}(mps.x1_idx) + sign(h_sl)*dh_sl;

% wrt p
Jfp_in = Jfp_mj(pi_ej(iq,1),x1,xd_ti{iq});
Jfp_out = Jfp_mj(pi_ej(iq,2),x1,xd_ti{iq});
dhp_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jfp_in+Jfp_out))...
    -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*(dh*(Jfp_in-Jfp_out));
Jsl_cond{iq}(mps.p_idx) =  Jsl_cond{iq}(mps.p_idx) + sign(h_sl)*dhp_sl;

for k = 1:nt % span all delays

        ik = id_ti{iq,k}(1); % index of containing segment
        Wdk = Wd_ti{iq,k}{1}; % interpolation coefficients
        WdTk = reshape(sum(WT_ti{iq,k}(1,:,:,:),4),dim,N); % derivatives wrt T
        WdPk = reshape(WP_ti{iq,k}(1,:,:),dim,Np); % derivatives wrt p
        mps_k = data.seg{ik}.coll_seg.maps;
        
        Jfk_in = Jf_mj(pi_ej(iq,1),x1,xd_ti{iq},k);  % incoming Jacobian matrix
        Jfk_out = Jf_mj(pi_ej(iq,2),x1,xd_ti{iq},k); % outgoing Jacobian matrix
        dhd_sl = 1/(dh*f_in - dh*f_out)*(dh*(Jfk_in+Jfk_out))...
            -(dh*f_in + dh*f_out)/(dh*f_in - dh*f_out)^2*...
            (dh*(Jfk_in-Jfk_out));

        % wrt xbp
        Jsl_cond{ik}(mps_k.xbp_idx) = Jsl_cond{ik}(mps_k.xbp_idx) ...
            + sign(h_sl)*dhd_sl*Wdk;

        % wrt p
        Jsl_cond{ik}(mps_k.p_idx) = Jsl_cond{ik}(mps_k.p_idx) ...
            + sign(h_sl)*dhd_sl*WdPk;

        % wrt xbp (for state dependent delays)
        if sys.sd_delay
            Wdxk = reshape(Wx_ti{iq,k}(1,:,:),dim,dim); % derivatives of tau wrt x
            Jsl_cond{iq}(mps.x1_idx) =  Jsl_cond{iq}(mps.x1_idx) ...
                + sign(h_sl)*dhd_sl*Wdxk;
         end
        
        % wrt T1
        for n = 1:N
            mps_n = data.seg{n}.coll_seg.maps;
            Jsl_cond{n}(mps_n.T_idx) = Jsl_cond{n}(mps_n.T_idx) ...
                + sign(h_sl)*dhd_sl*WdTk(:,n);
        end

        % Extra terms for neutral delayed terms
        if tau_type(k) == 2
            Jsl_cond{ik}(mps{ik}.T_idx) = Jsl_cond{ik}(mps{ik}.T_idx)...
                - sign(h_sl)*dhd_sl*Wdk*xbp{ik}/T(ik);
        end

end

% Return COCO compatible outputs
d = data;
f = abs(h_sl) - 1;
J = cell2mat(Jsl_cond);

end




