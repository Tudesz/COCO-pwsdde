function X = struct2func_tau(sys,p,ind,type)
%STRUCT2FUNC_TAU Convert time delay definitions from old struct data to new
%function definition form
% Input:
%   sys: old system definition data structure
%       -> tau: system time delays
%       -> tau_dp: derivative of time delays wrt system parameters
%   p: system parameter vector
%   ind: index of requested time delay
%   type: flag for which property is requested:
%       1, type of the delayed term:
%           -> 1, fix point delay (normal)
%           -> 2, fix point delay (neutral)
%       2, value of the time delay tau(p)
%       3, parameter jacobian of the time delay Jp_tau(p)
% Output:
%   X: appropriate evaluation of the delay or its parameter Jacobain

% Select vector field mode

% Select the requested time delay
switch type
    case 1 % Delay type
        X = 1;
    case 2 % Time delays
        tau = sys.tau(p);
        X = tau(ind);
    case 3 % Time delay derivatives
        tau_dp = sys.tau_dp(p);
        X = tau_dp(ind,:);
end

end

