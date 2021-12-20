% Romberg�㷨��Ͽ�����ų���
clear;
format long;
f = @(x) cosh(x/5); % ���ŵ����ߺ�����˫������
a = 0; 
b = 3;
tol = 1e-6;

h = b - a;
T1 = (h / 2) * (1 + f(b));
k = 1;
T2 = 0;
fprintf('T%d= %f\n', 2^0, T1);
flag = true;
while (flag)
    S = 0;
    x = a + h / 2;
    while (x < b)
        S = S + f(x);
        x = x + h;
    end
    T2 = T1 / 2 + h * S / 2;   % ��������������һ�����ֽ��
    fprintf('T%d= %f\n', 2^k, T2);
    S2 = T2 + (T2 - T1) / 3;  % �������������Simpson����
    fprintf('S%d= %f\n', 2^(k-1), S2);
    if k == 1
        k = k + 1;
        h = h / 2;
        T1 = T2;
        S1 = S2;    
    else
        C2 = S2 + (S2 - S1) / 15;   % �������������Cotes����
        fprintf('C%d= %f\n', 2^(k-2), C2);
        if k == 2
            C1 = C2;
            k = k + 1;
            h = h / 2;
            T1 = T2;
            S1 = S2;    
        else
            R2 = C2 + (C2 - C1) / 63;  % �������������Romberg����
            fprintf('R%d= %f\n', 2^(k-3), R2);
            if k == 3
                R1 = R2;
                C1 = C2;
                k = k + 1;
                h = h / 2;
                T1 = T2;
                S1 = S2;                
            else if abs(R2 - R1) >= tol
                    R1 = R2;
                    C1 = C2;
                    k = k + 1;
                    h = h / 2;
                    T1 = T2;
                    S1 = S2;
                else
                    flag = false;
                    fprintf('Romberg�㷨��Ͽ�����ų��� ��������%f��\n',R2 * 2);
                end
            end
        end
    end
end
