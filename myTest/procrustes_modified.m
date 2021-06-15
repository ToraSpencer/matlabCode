function [d, Z, transform] = procrustes_modified(point11, point2, varargin)
%,'scaling',0,'reflection',0

pnames = {   'scaling'  'reflection'};
dflts =  {       true         'best'};
 
[n, m]   = size(point11);
[ny, my] = size(point2);

% Center at the origin.
muX = mean(point11,1);
muY = mean(point2,1);
X0 = point11 - repmat(muX, n, 1);
Y0 = point2 - repmat(muY, n, 1);

ssqX = sum(X0.^2,1);
ssqY = sum(Y0.^2,1);
constX = all(ssqX <= abs(eps(class(point11))*n*muX).^2);
constY = all(ssqY <= abs(eps(class(point11))*n*muY).^2);
ssqX = sum(ssqX);
ssqY = sum(ssqY);

if ~constX && ~constY
 
    normX = sqrt(ssqX); % == sqrt(trace(X0*X0'))
    normY = sqrt(ssqY); % == sqrt(trace(Y0*Y0'))

    % Scale to equal (unit) norm.
    X0 = X0 / normX;
    Y0 = Y0 / normY;

    % Make sure they're in the same dimension space.
    if my < m
        Y0 = [Y0 zeros(n, m-my)];
    end

    % The optimum rotation matrix of Y.
    A = X0' * Y0;
    [L, D, M] = svd(A);
    T = M * L';
 
    haveReflection = (det(T) < 0);

    traceTA = sum(diag(D));  
    
    b = 1;
    d = 1 + ssqY/ssqX - 2*traceTA*normY/normX;

    if nargout > 1
        Z = normY*Y0 * T + repmat(muX, n, 1);
    end
 
    if nargout > 2
        if my < m
            T = T(1:my,:);
        end
        c = muX - b*muY*T;
        transform = struct('T',T, 'b',b, 'c',repmat(c, n, 1));
    end

elseif constX
    d = 0;
    Z = repmat(muX, n, 1);
    T = eye(my,m);
    transform = struct('T',T, 'b',0, 'c',Z);
    
else 
    d = 1;
    Z = repmat(muX, n, 1);
    T = eye(my,m);
    transform = struct('T',T, 'b',0, 'c',Z);
end