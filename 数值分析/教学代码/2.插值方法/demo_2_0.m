% ����ʽ��ֵ �����λ���ݵ�
% ��return����������
clear;clc;
xl = -3; xr = 3; yb= -3; yt = 3;    % ��������
plot([xl xr], [0 0], 'k', [0 0], [yb yt], 'k');
grid on;
xlist = []; ylist = [];    % �洢�û�ѡ���ֵ�ڵ�ĺ�������
k = 0;
while(0 == 0)
    [xnew, ynew] = ginput(1);   % �û������������ѡ���ֵ�ڵ�
    if length(xnew) < 1
        break
    end
    k = k+1;
    xlist(k) = xnew;
    ylist(k) = ynew;
    [m,c] = d_d(xlist, ylist);   % ���ݲ�ֵ�ڵ������̱����ضԽ���Ԫ��
    
    poly_exp = 'y = '; omiga_exp = '';
    for index = 1 : length(c)  
        if c(index) >= 0
            coef_str = ['+', num2str(c(index))];
        else
            coef_str = num2str(c(index));
        end
      poly_exp = [poly_exp, coef_str,omiga_exp];
        if -xlist(index) >= 0 
            omiga_sub = ['+', num2str(-xlist(index))];
        else
            omiga_sub = num2str(-xlist(index));
        end
      omiga_exp = [omiga_exp, '(x', omiga_sub,')'];
    end
    
    
    x = xl : 0.01 : xr;          % ����ֵ�㣬�������ߵĺ�����
    y = nest(c, x, xlist);       % �������ֵ���������    
    plot(xlist, ylist, 'o', x, y, [xl xr], [0 0], 'k', [0 0], [yb yt], 'k');
    title('����ʽ��ֵ��ʾ');
    text(xl, yt - 0.5, poly_exp);
    axis([xl xr yb yt]); grid on;
end
    