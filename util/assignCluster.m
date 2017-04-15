function Z=assignCluster(Y,X)
%
% Assining the same clasters as the nearest initial points 
% to irregular points.
%  

NaNind=find(isnan(sum(Y,1))==1); NaNdata=X(:,NaNind);

dist=repmat(sum(X.^2,1),[length(NaNind),1])...
    +repmat(sum(NaNdata.^2,1)',[1,size(X,2)])...
    -2*(NaNdata')*X;
dist(:,NaNind)=NaN; % Removing distance to itself.

[~,minind]=min(dist,[],2);

Z=Y; Z(:,NaNind)=Y(:,minind);