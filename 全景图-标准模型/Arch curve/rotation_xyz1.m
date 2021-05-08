function [ixyz] = rotation_xyz1(ixyz,thx,thy,thz,pt,ang_type)
%%
%% rotate_xyz? ?in angle thx (along x-axis)
%%? ?? ?? ?? ???in angle thy (along y-axis)
%%? ?? ?? ?? ???in angle thz (along z-axis)
%% pt? ?? ?? ???center of rotation
%% ang_type? ???type of the angle. 1: degree, else: rad? ?? ?? ?? ?


if ang_type == 1
    thx = thx*pi/180;
    thy = thy*pi/180;
    thz = thz*pi/180;
end

ixyz(:,1) = ixyz(:,1) - pt(1);
ixyz(:,2) = ixyz(:,2) - pt(2);
ixyz(:,3) = ixyz(:,3) - pt(3);

mz = [ cos(thz) -sin(thz) 0
    sin(thz) cos(thz) 0
    0  0 1 ];
my = [ cos(thy) 0 sin(thy) 
    0  1 0 
 -sin(thy) 0 cos(thy) ];
mx = [ 1 0  0
    0 cos(thx) -sin(thx) 
    0 sin(thx) cos(thx) ];
ixyz = ixyz*mx;
ixyz = ixyz*my;
ixyz = ixyz*mz;

ixyz(:,1) = ixyz(:,1) + pt(1);
ixyz(:,2) = ixyz(:,2) + pt(2);
ixyz(:,3) = ixyz(:,3) + pt(3);