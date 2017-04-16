clearvars; close all; clc;
addpath('util');
addpath('Manifolds');

rand_id     = 1;
nc          = 50;
n           = nc*3;
dim1        = 5;
dim2        = 2;

settings_id_list = [1, 4];

rng(rand_id)

%% set options
op.samples=n; % Num. of samples.
op.bnum=min(100, op.samples); % number of basis functions
op.cvfold=5; % N-fold cross validation
op.sigma_list=logspace(-3,1,8);
op.lambda_list=logspace(-3,1,8);
op.maxiter=200; % maximum number of iteration in clustering update
op.tol=1e-4; % Stopping criteria in clustering
op.basis='DerGaussian';
op.regularizer='L2';
op.updateform='vector';
op.manifolds='Grassmann';
r = 1e-3;

figure(1); clf;
for ite = 1:length(settings_id_list)
    settings_id = settings_id_list(ite);
    [A,X,true_clusters,op] = toydata(rand_id,settings_id,op,dim1,dim2);

    method = 'R-LSLDGC';
    fprintf('R-LSLDGC\n');
    tic_id = tic;
    [sigma, lambda, C, AC_dist, P, l] = CV_RLSLDG(A, op);
    theta_RLSLDG = computeTheta(sigma, lambda, AC_dist, P, l, op);
    Y.RLSLDGC=RLSLDGClust(A, sigma, theta_RLSLDG, C, op); 
    time = toc(tic_id);

    % evaluation
    [Nclusters,~,ll]=econncomp_geo(Y.RLSLDGC,r);
    ARI = valid_RandIndex(ll, true_clusters);

    geodist = zeros(op.samples,op.samples);
    for i = 1 : op.samples
        for j = 1:op.samples
            geodist(i,j) = calc_squared_geodist(A(:, :, i), A(:, :, j), op);
        end
    end
    subplot(2, 3, 1 + (ite-1)*3);
    imagesc(geodist); %colorbar;
    axis equal tight
    title('Riemannian distance');
    if ite == 1
        ylabel('\gamma=0');
    else
        ylabel('\gamma=\pi/2');
    end

    eucdist = zeros(op.samples,op.samples);
    for i = 1 : op.samples
        for j = 1:op.samples
            eucdist(i,j) = sum(sum((X(:,i)-X(:,j)).^2,1),2);
        end
    end
    subplot(2, 3, 2 + (ite-1)*3);
    imagesc(eucdist); %colorbar;
    axis equal tight
    title('Euclidean distance');

    ret_geodist = zeros(op.samples,op.samples);
    for i = 1 : op.samples
        for j = 1:op.samples        
            ret_geodist(i,j) = calc_squared_geodist( ...
                Y.RLSLDGC(:, :, i), Y.RLSLDGC(:, :, j), op);
%             op.dim2 - trace(Y.RLSLDGC(:,:,i)'*Y.lsldggeo(:,:,j)*Y.lsldggeo(:,:,j)'*Y.lsldggeo(:,:,i));
        end
    end
    subplot(2, 3, 3 + (ite-1)*3);
    imagesc(ret_geodist); %colorbar;
    axis equal tight
    title('Result of R-LSLDGC');
    xlabel(sprintf('ARI: %.2f\n', ARI));
    
    fprintf('time: %.2f [sec.]\n', time);
end

print('-dpng', 'result.png', '-r0');

