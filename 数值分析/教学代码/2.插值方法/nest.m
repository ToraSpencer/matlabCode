% ţ�ٲ�ֵ����ʽ��ֵ
% ������� c ��ţ�ٲ�ֵ����ʽ��ϵ����Ҳ���ǲ��̱�ĶԽ���
% ������� x ����ֵ���Ա���ȡֵ������
% ������� xlist ��֪�Ĳ�ֵ�ڵ�
% ������� y ��ֵ�������
function y = nest(c, x, xlist)
n = length(c);
result = 0;
item = 1;
for k = 1 : n  
    result = result + c(k) * item;
    item = item .* (x - xlist(k));   % �����ĵ��
end
y = result;