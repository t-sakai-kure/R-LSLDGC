function [A,X,true_clusters,op] = toydata(rand_id,settings_id,op,dim1,dim2)
% A : samples (dim1 by dim2 by samples)
% X : vectorize samples A  (dim1*dim2 by samples)

rng(rand_id);
op.dim1 = dim1;
op.dim2 = dim2;
op.dim = op.dim1*op.dim2;
op.bnum = min(100,op.samples);

D = randn([op.dim1,op.dim2]);    
[U,~,~] = svd(D);
V = U(:,1:op.dim2);

%%

R1 = zeros(op.dim2,op.dim2,op.samples);
sigma = [0 pi/8 pi/4 pi/2];
th = sigma(settings_id)*randn(op.samples, 1);

if op.dim2 == 2
    for i = 1: op.samples
        R1(:,:,i) =  [cos(th(i)) -sin(th(i));
            sin(th(i)) cos(th(i))] ;
    end
elseif op.dim2 == 3
    for i = 1: op.samples
        R1(:,:,i) =  [cos(th(i)) -sin(th(i)) 0;
            sin(th(i)) cos(th(i)) 0;
            0 0 1] ;
    end
end

t = [pi/15*randn(1/3*op.samples,1); 2*pi/3 + pi/15*randn(1/3*op.samples,1); 4*pi/3 + pi/15*randn(1/3*op.samples,1)];
ROT = zeros(2,2,op.samples);
for i = 1:op.samples
    ROT(:,:,i) = [cos(t(i)) -sin(t(i));
        sin(t(i)) cos(t(i))];
end

A = zeros(op.dim1,op.dim2, op.samples);
if op.dim1 == 2
    for i = 1 : op.samples
        A(:,:,i) = ROT(:,:,i)*V(:,:,1)*R1(:,:,i);
    end
else
    for i = 1: op.samples
        A(:,:,i) =  blkdiag(ROT(:,:,i),eye(op.dim1-2)) *V(:,:,1)*R1(:,:,i);
    end
end

%%
X = reshape(A,[op.dim,op.samples]);
true_clusters = [ones(op.samples/3, 1); 2*ones(op.samples/3, 1); 3*ones(op.samples/3,1)];
end

