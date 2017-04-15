% geod = geoddist(X[,Y,w]) Matrix of squared (weighted) geodesic distances d(X,Y)


function geod = geodist(X,Y)
% for Grassmann manifolds

if nargin==1	% Fast version for common case sqdist(X) 
    [D1,D2] = size(X);
    tX = permute(X,[2,1,3]);
    geod = D2-trace(tX*X*tX*X);
  return
end

% ---------- Argument defaults ----------
%if ~exist('Y','var') | isempty(Y) Y = X; eqXY = 1; else eqXY=0; end;

% The intervector squared distance is computed as (x-y)² = x²+y²-2xy.
% We ensure that no value is negative (which can happen due to precision loss
% when two vectors are very close).
[D1,D2] = size(X);

tX = permute(X,[2,1,3]);
tY = permute(Y,[2,1,3]);

geod = zeros(size(X,3),size(Y,3));
for i = 1:size(X,3)
    for j = 1:size(Y,3)
         geod(i,j) = D2-trace(tX(:,:,i)*Y(:,:,j)*tY(:,:,j)*X(:,:,i));        
    end
end


