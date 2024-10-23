function [data, sol] = pwsdde_coll_init(prob,src_data,NTST)
%PWSDDE_COLL_INIT Initialize a collocation segment using the 'coll'
%COCO toolbox (based on 'ode_init_data', 'coll_add', and
%'coll_construct_seg')
% Input: 
%   prob: COCO problem structure to append new functions to
%   src_data: initialization source data, must contain the fields
%    -> t0: initiali time mesh (N x 1 column vector)
%    -> x0: solution guess on t0 (N x dim matrix)
%    -> p0: default parameter vector
%    -> pnames: cell array of parameter names
%   NTST: starting mesh segment number (default 10)
% Output:
%   data: data on the mesh used by the COCO 'coll' toolbox
%   sol: solution data structure defined on the new mesh

% load default options
data = coco_func_data('protect');
data.oid = 'seg';
tbid = coco_get_id(data.oid, 'coll');
data.coll = coll_get_settings(prob, tbid, []);

% append the fields of src_data
data.t0 = src_data.t0;
data.x0 = src_data.x0;
data.p0 = src_data.p0;
data.pnames = src_data.pnames;

% adjust orientation of x0
if size(data.x0,1) ~= length(data.t0)
    data.x0 = data.x0.';
end

 % initialize solution structure
[sol, data] = coll_read_solution('', '', data);
if nargin>2
    data.coll.NTST = NTST;
end
[data, sol] = coll_init_data(data, sol);
data = data.data;
data = rmfield(data,'oid');

end

