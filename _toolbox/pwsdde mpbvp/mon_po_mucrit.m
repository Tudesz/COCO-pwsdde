function [d, f] = mon_po_mucrit(prob,data,u)
%MON_PO_MUCRIT Monitor the critical Floquet multiplier of the periodic 
% orbits
% Input:
%   prob: COCO problem structure
%   data: collocation data structure
%   u: vector of continuation variables
% Output: 
%   d: collocation data structure (copy of the input "data")
%   f: critical Floquet multiplier

% Formulate the monodromy matrix
Th = pwsdde_mon(data,u);

% Evaluate critical Floquet-multiplier from monodormy matrix
mu = eig(full(Th));
[~,ind] = max(abs(mu));
muc = mu(ind);

% COCO compatible outputs
d = data;
f = [muc; abs(muc)];

end

