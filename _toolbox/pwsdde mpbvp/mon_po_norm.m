function [d, f] = mon_po_norm(prob,data,u)
%MON_PO_AMPL Monitor the norm of the vector of state variables as well as
%the amplitudes of the state variables
% Input:
%   prob: COCO problem structure
%   data: collocation data structure
%   u: vector of continuation variables
% Output: 
%   d: collocation data structure (copy of the input "data")
%   f: evaluation of norm((max(u_i)-min(u_i))/2)

% Unpack solution and segment meshes
N = length(data.seg);
xbp = cell(1,N);    % state at mesh nodes
for i=1:N
    mps = data.seg{i}.coll_seg.maps;
    uidx = data.seg{i}.uidx;
    xbp{i} = reshape(u(uidx(mps.xbp_idx)),mps.xbp_shp);
end

% Calculate vibration amplitudes
x = cell2mat(xbp);
ampl = (max(x,[],2)-min(x,[],2))./2;

% COCO compatible outputs
d = data;
f = [norm(x); ampl];

end