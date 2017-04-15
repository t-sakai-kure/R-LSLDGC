function [sigma, lambda, C, AC_dist2, P, l] = CV_RLSLDG(A, op, C)
%
% Estimating a Log-Density Gradient on Reimaniann manifolds
%
% X: (sample) by (dim) matrix
% op: options
% C: center points of kernelsf
%

narginchk(2,3);
                                                                                                               
if nargin < 3 || isempty(C)
    cind = randperm(op.samples,op.bnum);
    C = A(:,:,cind);
else
    cind = [];
end

% Distance from centers (size: bnum by samples)
AC_dist2 = zeros([op.bnum op.samples]);
L       = cell(op.bnum, op.samples);
for k = 1:op.bnum
    for i= 1:op.samples
        AC_dist2(k,i) = calc_squared_geodist(A(:, :, i), C(:, :, k), op);
        L{k, i} = calc_logmap(A(:, :, i), C(:, :, k), op);
    end
end

P = cell(op.dim, 1);
l = cell(op.dim, 1);
for dd = 1:op.dim
    P{dd} = zeros(op.bnum, op.samples);
    l{dd} = zeros(op.bnum, op.samples);
end


for k = 1:op.bnum
    fprintf('bnum: %d/%d\n', k, op.bnum);
    M1 = C(:,:,k)*C(:,:,k)';
    for i = 1:op.samples        
        M2 = A(:,:,i)'*M1*A(:,:,i);
        M3 = M1*A(:, :, i);
        M4 = A(:,:,i)*A(:,:,i)'*M1;
        
        for dd = 1:op.dim            
            Fx = zeros(op.dim1, op.dim2);
            for j = 1:op.dim
                [m,  n]  = ind2sub([op.dim1 op.dim2], dd);
                [m2, n2] = ind2sub([op.dim1 op.dim2], j);
                Fx(j) = M1(m, m2)*(n == n2) - M2(n, n2)*(m == m2) ...
                    - A(m, n2, i)*M3(m2, n) - M4(m, m2)*(n == n2);
            end
            rgrad = (eye(op.dim1) - A(:, :, i)*A(:, :, i)')*Fx;
            P{dd}(k, i) = rgrad(dd);
            l{dd}(k, i) = L{k, i}(dd);
        end
    end
end
clear L

% cross validation
cv_fold  = (1:op.cvfold);
cv_split = floor((0:op.samples-1)*op.cvfold./op.samples)+1;
cv_index = cv_split(randperm(op.samples));

%%
score_cv = zeros(length(op.sigma_list), length(op.lambda_list), length(cv_fold));
for sigma_index = 1:length(op.sigma_list)
    sigma_tmp = op.sigma_list(sigma_index);
    GauKer = exp(-AC_dist2/(2*sigma_tmp^2)); % bnum by samples
           
    for kk = cv_fold
        G_tr = zeros(op.bnum, op.bnum);
        G_te = zeros(op.bnum, op.bnum);
        h_tr = zeros(op.bnum, 1);
        h_te = zeros(op.bnum, 1);        
        for dd = 1:op.dim            
            phi = P{dd}.*GauKer + (l{dd}.^2).*GauKer/sigma_tmp^2;
            psi_tr  = l{dd}(:, cv_index~=kk).*GauKer(:, cv_index~=kk);
            psi_te  = l{dd}(:, cv_index==kk).*GauKer(:, cv_index==kk);
            
            phi_tr  = phi(:, cv_index~=kk);   %bnum by index~=kk
            phi_te  = phi(:, cv_index==kk);   %bnum by index==kk
            
            G_tr = G_tr + psi_tr*psi_tr'/size(psi_tr,2);
            G_te = G_te + psi_te*psi_te'/size(psi_te,2);
            h_tr = h_tr + mean(phi_tr, 2);
            h_te = h_te + mean(phi_te, 2);
        end % dd
        for lambda_index = 1:length(op.lambda_list)
            lambda_tmp = op.lambda_list(lambda_index);
            thetah = -(G_tr + lambda_tmp*eye(op.bnum))\h_tr;
            term = thetah'*G_te*thetah + 2*thetah'*h_te;
            score_cv(sigma_index, lambda_index, kk) = term;
        end % lambda
    end % cv
end % sigma

[score_cv_tmp, lambda_index] = min(mean(score_cv, 3), [], 2);

if all(score_cv_tmp(:)==Inf)
     disp('ALL Inf')
end
[~, sigma_index] = min(score_cv_tmp);
lambda = op.lambda_list(lambda_index(sigma_index));
sigma  = op.sigma_list(sigma_index);
fprintf('sigma=%g, lambda=%g\n', sigma, lambda);

end
