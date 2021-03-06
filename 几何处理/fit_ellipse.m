% 最小二乘法拟合椭圆
% Input: x0,y0 are lists of coordinates [x0, y0]是所有样本点；
function x = fit_ellipse(x0, y0,type)
    if nargin == 2
        type = 'general';
    end
    
    switch type
        case 'general'  % 一般椭圆方程a*x^2 + b*x*y + c*y^2+d*x+e*y+f = 0;
            % x = [a, b, c, d, e, f]，椭圆方程写为A*x = 0;
            A = [ x0.*x0, x0.*y0, y0.*y0, x0, y0, ones(size(x0)) ];
            
            % Build scatter matrix
            S = A'*A;
            
            % Build 6x6 constraint matrix
            B(6,6) = 0; 
            B(1,3) = 2; 
            B(2,2) = -1; 
            B(3,1) = 2;
            
            % Solve eigensystem
            [gevec, geval] = eig(inv(S)*B);
            
            % Find the positive eigenvalue
            [PosR, PosC] = find(geval > 0 & ~isinf(geval));
            
            % Extract eigenvector corresponding to positive eigenvalue
            x = gevec(:,PosC);
            
        case 'standard'            % 标准椭圆方程a*x^2 + c*y^2+d*x+e*y+f = 0; 没有旋转；
            % x = [a, c, d, e, f]，椭圆方程写为A*x = 0; 约束条件：a*c > 0，令a*c = 4 
            A = [ x0.*x0, y0.*y0, x0, y0, ones(size(x0)) ];
            S = A'*A;
            B(5,5) = 0; 
            B(1,2) = 2; 
            B(2,1) = 2;
            [gevec, geval] = eig(inv(S)*B);
            [PosR, PosC] = find(geval > 0 & ~isinf(geval));
            x = gevec(:,PosC);
            
        otherwise
            error('type error!');
    end
end
