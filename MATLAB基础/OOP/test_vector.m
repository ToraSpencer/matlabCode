clc;
clear all;
%%
v0 = vector2D();
disp(v0);
fprintf('rou = %g\n',v0.rou);

%%
v1 = vector2D(1,2);
disp(v1);


%%
v1.x = 100;
disp(v1);

%%
v1.normalize();
disp(v1);


%%
v2 = vector3D(1,2,3);
disp(v2);

%%
v2.normalize();
disp(v2);

