function [W, Wp, Wpp] = pwsdde_coll_W(maps,tq)
%PWSDDE_COLL_W Interpolation coefficients for a given point considering
%the discretization of COCO (using 'coll_L', 'coll_Lp' and 'coll_Lpp')
%   maps: coll_seg maps data structure
%   tq: vector of querry points rescaled to the [0 1] interval (in column
%   vector form)
% Output:
%   W: Matrix of Lagrange interpolation coefficients
%   Wp: Lagrange interpolation coefficients of the first derivatives
%   Wpp: Lagrange interpolation coefficients of the second derivatives

% initialization
NTST = maps.NTST;
NW = size(maps.W,2);
dim = maps.x_shp(1);
NCOL = maps.x_shp(2)/NTST;

% find contiaing mesh subintervals
tq_N = ceil(NTST*tq);
tq_N(tq_N==0) = 1;
tq_N(tq_N>NTST) = NTST;
tq_dN = reshape(repmat(tq_N,1,dim).',[],1);

% fill up interpolation matricies
W = zeros(length(tq_dN),NW);    % xq wrt xbp
Wp = zeros(length(tq_dN),NW);   % dxq wrt xbp
Wpp = zeros(length(tq_dN),NW);  % ddxq wrt xbp
tm = linspace(-1, 1, NCOL+1).';
for i = reshape(unique(tq_N),1,[])
    tr = tq(tq_N==i)*NTST*2 - 2*i + 1; % NTST*2 multiplier!
    pmap = coll_L(tm, tr);
    dmap = coll_Lp(tm, tr);
    ddmap = coll_Lpp(tm, tr);
    col_idx = (i-1)*dim*(NCOL+1)+1:i*dim*(NCOL+1);
    W(tq_dN==i,col_idx) = kron(pmap, eye(dim));
    Wp(tq_dN==i,col_idx) = kron(dmap, eye(dim));
    Wpp(tq_dN==i,col_idx) = kron(ddmap, eye(dim));
end

% scale with the multiplier of tr
Wp = 2*NTST*Wp;
Wpp = (2*NTST)^2*Wpp;

end



