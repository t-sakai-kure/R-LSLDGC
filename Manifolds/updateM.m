function newZ = updateM(Z, M, op)

manifolds = lower(op.manifolds);

switch manifolds
    case 'grassmann'
        [U,S,V] = svd(M,0);
        cosS = diag(cos(diag(S)));
        sinS = diag(sin(diag(S)));
        
        Z = (Z*V*cosS + U*sinS) *V';
        [newZ, unused] = qr(Z, 0); %#ok
        
%     case 'lie'
%          newZ = Z*expm(Z\M);
         
%     case 'stiefel'
%         [Q,R] = qr((eye(op.dim1)-Z*Z')*M);
%         BD = expm([E -R'; R zeros(op.dim2,op.dim2)])*[eye(op.dim2);zeros(op.dim2,op.dim2)];
%         B = BD(1:op.dim2,:);
%         D = BD(op.dim2+1:2*op.dim2,:);
%         newZ = Z * B + Q * D;
        
%     case 'spd'
%         newZ = symm(Z*real(expm(Z\(1.0*M))));
        
    otherwise
        error('%s is not supported', manifolds);
        
        
end

end