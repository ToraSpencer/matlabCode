% ����������ֵ �����λ��ֵ��
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
    if k == 1
        continue;
    end
    x = xl : 0.01 : xr;          % ����ֵ�㣬�������ߵĺ�����   
    pp = spline(xlist, [0, ylist, 0]);  % ��һ��߽�����
    %pp = spline(xlist, ylist);  % not a knot��������
    y = ppval(pp, x);                   % �������ֵ���������    
    plot(xlist, ylist, 'o', x, y, [xl xr], [0 0], 'k', [0 0], [yb yt], 'k');   
    title('����������ֵ��ʾ');
    axis([xl xr yb yt]); grid on;
end
    