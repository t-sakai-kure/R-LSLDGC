function theta=computeTheta_geo(sigma,lambda,op,AC_dist, P, l)
%
% Computing theta
%
narginchk(6,6);

GauKer=exp(-AC_dist/(2*sigma^2));
G = zeros(op.bnum, op.bnum,op.dim);
h = zeros(op.bnum,op.dim);
for dd=1:op.dim
    fprintf('d: %d/%d\n', dd, op.dim)
    phi = P{dd}.*GauKer + (l{dd}.^2).*GauKer/sigma^2;
    psi=l{dd}.*GauKer; %bnum *samples
    G(:,:,dd)=psi*psi'/size(psi,2); % bnum * bnum
    h(:,dd)=mean(phi,2);
end

sumG = sum(G,3);
sumh = sum(h,2);
theta = -(sumG+lambda*eye(op.bnum))\sumh;

theta=repmat(theta',[op.dim,1]); % dim * bnum

end