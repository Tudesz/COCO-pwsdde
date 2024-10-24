function sys = struct2func(struct)
%STRUCT2FUNC Convert old struct system definition to new function form
% Input:
%   struct: old system definition data structure
%     -> f{m}: governing autonomous vector spaces R^n -> R^n (DDE rhs)
%     -> h{i}: switching or event surfaces R^n -> R
%     -> g{i}: event maps R^n->R^n
%     -> Ju{m}: vector space Jacobians R^n -> R^(n x n)(DDE rhs)
%     -> Jud{m,k}: vector space Jacobians wrt x(t-tau(k))
%     -> dh{i}: event surface Jacobians R^n -> R^(1 x n)
%     -> dg{i}: event map Jacobians R^n -> R^(n x n)
%     -> dhd{i,k}: event surface Jacobians wrt x(t-tau(k))
%     -> dgd{i,k}: event map Jacobians wrt x(t-tau(k))
%     -> Jp{m}: Jacobian of f R^n -> R^(n x l)
%     -> dhp{i}: Jacobian of h R^n -> R^(1 x l)
%     -> dgp{i}: Jacobian of g R^n -> R^(n x l)
%     -> tau: system time delays
%     -> tau_dp: derivative of time delays wrt system parameters
%     -> pm{i}: mode changes on event surfaces
% Output:
%   sys: functions compatible with the new continuation framework
%     -> f: vector field definitions
%     -> e: event functions and maps
%     -> tau: time delays
%     -> mode_no: number of vector field modes
%     -> event_no: number of possible events

if isfield(struct,'tau')
    sys.tau_no = size(struct.Jud,2);
else
    sys.tau_no = 0;
end
sys.mode_no = length(struct.f);
sys.event_no = length(struct.h);

sys.f = @(x,xd,p,mode,type,l) struct2func_f(struct,x,xd,p,mode,type,l);
sys.e = @(x,xd,p,id,type,l) struct2func_e(struct,x,xd,p,id,type,l);
if sys.tau_no > 0
    sys.tau = @(p,ind,type) struct2func_tau(struct,p,ind,type);
else
    sys.tau = @(p,ind,type) [];
end



end


