%delete(gcp):�رյ�ǰ����ѭ����
%parpool(���õļ��������Ŀ)���򿪲���ѭ����

clc;
clear all;

delete(gcp);
K = 1000000;

%%
sum = 0;
tic;
for i = 1:K
    sum = sum+i;
end
t0 = toc;

%%
sum = 0;
parpool(1);
tic;
parfor i = 1:K
    sum = sum+i;
end
t1 = toc;

%%
sum = 0;
delete(gcp);
parpool(2);
tic;
parfor i = 1:K
    sum = sum+i;
end
t2 = toc;

%%
sum = 0;
delete(gcp);
parpool(3);
tic;
parfor i = 1:K
    sum = sum+i;
end
t3 = toc;

%%
sum = 0;
delete(gcp);
parpool(4);
tic;
parfor i = 1:K
    sum = sum+i;
end
t4 = toc;


