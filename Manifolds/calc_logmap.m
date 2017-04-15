function logmap = calc_logmap(A, C, op)

manifolds = lower(op.manifolds);

switch manifolds
    case 'grassmann'
        logmap = (eye(op.dim1) - A*A')*(C*C')*A;
        
%     case 'lie'
%         logmap = A*logm(A\C);
        
%     case 'stiefel'
%         logmap = A*C'*A - C;
        
%     case 'spd'
%         logmap = symm(A*real(logm(A\C)));
        
    otherwise
        error('%s is not supported', manifolds);
end

end