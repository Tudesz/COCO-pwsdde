function [data, sol] = coll_init_data(data, sol)
%COLL_INIT_DATA   Initialize the data structure of the 'coll' segment

t0 = sol.tbp-sol.T0;

NCOL = data.coll.NCOL;
seg.int  = coll_interval(NCOL, data.xdim);

NTST = data.coll.NTST;
seg.maps = coll_maps(seg.int, NTST, data.pdim);

if abs(sol.T)>eps
  t  = linspace(0, NTST, numel(t0));
  tt = interp1(t, t0, 0:NTST, 'linear');
  tt = tt*(NTST/tt(end));
else
  tt = 0:NTST;
end
seg.mesh = coll_mesh(seg.int, seg.maps, tt);

seg.fid  = coco_get_id(data.oid, 'coll');
data.coll_seg = seg;

x0 = sol.xbp;

if abs(sol.T)>eps
  t0 = t0/sol.T;
  [t0,uind] = unique(t0);
  x0 = x0(uind,:);
  x0 = interp1(t0, x0, seg.mesh.tbp)';
else
  x0 = repmat(x0(1,:), size(seg.mesh.tbp))';
end

sol.u0 = [x0(:); sol.T0; sol.T; sol.p];
  
if ~isempty(sol.t0)
  x0_t = sol.xbp_t0;
  if abs(sol.T)>eps
    x0_t = interp1(t0, x0_t, seg.mesh.tbp)';
  else
    x0_t = repmat(x0_t(1,:), size(seg.mesh.tbp))';
  end
  sol.t0 = [x0_t(:); sol.T0_t0; sol.T_t0; sol.p_t0];
end

end
