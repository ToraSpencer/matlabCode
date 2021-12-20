% ���Ƽ��� �㷨�ȶ��Ե���ʾ
clear;

% �ⷨһ   n��С����
S0 = 0.182;
S(1) = 1 - 5*S0;
for n = 2 : 8
    S(n) = 1/n - 5*S(n-1);
    S(n) = vpa(S(n),3);
end
disp('�ⷨһ��')
disp(S)

% �ⷨ��   n�Ӵ�С
S(8) = 0.0204;
for n = 8: -1 : 2
    S(n-1) = 1/(5*n) - S(n)/5;
    S(n-1) = vpa(S(n-1),3);
end
disp('�ⷨ����')
disp(S)

% ��ȷֵ
syms x
for n = 1 : 5 % 8
    St(n) = int(x^n/(x+5),0,1);
    St(n) = vpa(St(n),3);
end
disp('��ȷֵ��')
disp(S)