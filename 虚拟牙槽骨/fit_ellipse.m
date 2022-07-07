% Reference: Direct Least Squares Fitting of Ellipses
% Input: x,y are lists of coordinates
function a = fit_ellipse(x,y,type)
    if nargin == 2
        type = 'general';
    end
    switch type
        case 'general'
            % Build design matrix
            D = [ x.*x x.*y y.*y x y ones(size(x)) ];
            % Build scatter matrix
            S = D'*D;
            % Build 6x6 constraint matrix
            C(6,6) = 0; C(1,3) = 2; C(2,2) = -1; C(3,1) = 2;
            % Solve eigensystem
%             [gevec, geval] = eig(S,C);
            [gevec, geval] = eig(inv(S)*C);
            % Find the positive eigenvalue
            [PosR, PosC] = find(geval > 0 & ~isinf(geval));
            % Extract eigenvector corresponding to positive eigenvalue
            a = gevec(:,PosC);
        case 'standard'
            D = [ x.*x y.*y x y ones(size(x)) ];
            S = D'*D;
            C(5,5) = 0; C(1,2) = 2; C(2,1) = 2;
            [gevec, geval] = eig(inv(S)*C);
            [PosR, PosC] = find(geval > 0 & ~isinf(geval));
            a = gevec(:,PosC);
        otherwise
            error('type error!');
    end
end
