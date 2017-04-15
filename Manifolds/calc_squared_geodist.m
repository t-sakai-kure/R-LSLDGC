function geodist2 = calc_squared_geodist(A, C, op)

manifolds = lower(op.manifolds);

switch manifolds
    case 'grassmann'
        geodist2 = op.dim2 - trace(C'*(A*A')*C);
        
%     case 'lie'
%         geodist2 = sum(sum((logm(C\A)).^2,1),2);
        
%     case 'stiefel'
%         geodist2 = op.dim2-trace(A'*C);

%     case 'spd'
%         geodist2 = real(trace((real(logm(A\C)))^2));
    otherwise
        error('%s is not supported', manifolds);
end

end

