function A = coll_Lpp(ts, tz)
% Second derivatives of Lagrange interpolation coefficients

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1 1]), [1 q q q q]);
sj = repmat(reshape(ts, [1 q 1 1 1]), [p 1 q q q]);
sk = repmat(reshape(ts, [1 1 q 1 1]), [p q 1 q q]);
sl = repmat(reshape(ts, [1 1 1 q 1]), [p q q 1 q]);
sh = repmat(reshape(ts, [1 1 1 1 q]), [p q q q 1]);

t3 = sj(:,:,:,1,1) - sk(:,:,:,1,1);
t4 = sj(:,:,:,:,1) - sl(:,:,:,:,1);
t5 = zi - sh;
t6 = sj - sh;

idx1 = find(abs(t6)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(t4)<=eps);
idx4 = find(abs(sl(:,:,:,:,1)-sk(:,:,:,:,1))<=eps);

idx5 = find(abs(sl-sh)<=eps);
idx6 = find(abs(sk-sh)<=eps);

t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

t4(union(idx3,idx4)) = 1;
t4                   = 1.0./t4;
t4(union(idx3,idx4)) = 0;

t6(union(idx1, union(idx5, idx6))) = 1;
t5(union(idx1, union(idx5, idx6))) = 1;

A = sum(t3.*sum(t4.*prod(t5./t6,5),4),3);
end