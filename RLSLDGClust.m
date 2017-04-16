function Y=RLSLDGClust(A, sigma, theta, C, op)
%
% Clustering via mode seeking based on the LSLDG estimator.
%
narginchk(5,5);

theta = reshape(theta,[op.dim1,op.dim2,op.bnum]);
Y = zeros(size(A));

for ii=1:op.samples
    Z=A(:,:,ii); % fix one point

    iter=0; criteria=1; NaNcall=0;
    while criteria
        iter=iter+1;
        ZC_geodist = zeros(1, op.bnum);   %1 by bnum
        for i=1:op.bnum
            ZC_geodist(i) = calc_squared_geodist(Z, C(:, :, i), op);            
        end
        
        GauKer = repmat(exp(-ZC_geodist/(2*sigma^2)), [op.dim 1]); % dim x bnum
        GauKer = reshape(GauKer,[op.dim1,op.dim2,op.bnum]);
        W = zeros(op.dim1, op.dim2, op.bnum);
        for k = 1:op.bnum
            W(:,:,k) = calc_logmap(Z, C(:, :, k), op);
        end
        M= sum(theta.*W.*GauKer,3)./sum(theta.*GauKer,3);
       
        % convergence criteria
        d=sqrt(sum(sum(M.^2,1),2));
        % d=max(d(d<inf),[],2);
        criteria = (d > op.tol) & (iter < op.maxiter);
        fprintf('ii=%g, iter=%g, d=%g\n',ii,iter,d);
        
        Z = updateM(Z, M, op);
        
        % convergence criteria
    %    d=sqrt(op.dim2 - trace(Zold'*Z*Z'*Zold));
    %    criteria = (d > op.tol) & (iter < op.maxiter);
    %    fprintf('ii=%g, iter=%g, d=%g \n',ii,iter,d);
        if sum(isnan(sum(sum(Z,1))))~=0 || sum(isinf(sum(sum(Z,1))))~=0
            NaNcall=1;
            break;
        end
    end

    if NaNcall;
        Y(:,:,ii)=NaN;
    else
        Y(:,:,ii)=Z;
    end
    fprintf('\n');
end


infind=find(isinf(sum(sum(abs(Y),1),2))==1);
if ~isempty(infind)
    Y(:,infind)=NaN;
end

%fprintf('NaNpoints=%g\n',sum(isnan(sum(Y,1))));

% Assigning clusters to NaN points
if sum(isnan(sum(sum(Y,1),2)))~=0
    Y=assignCluster_geo(Y,A);
end

end

