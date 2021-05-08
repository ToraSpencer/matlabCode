function Pr=rot3d1(P,origin,dirct,theta)
% 将坐标点P绕着，过origin点，方向为dirct的直线，旋转theta角
% P：需要旋转的做标集合，n×3矩阵
% origin：转轴通过的点，1×3向量
% direct：转轴方向向量，1×3向量
% theta：旋转角度，单位弧度


dirct=dirct(:)/norm(dirct);
T = [dirct(1).^2*(1-cos(theta))+cos(theta)  dirct(1)*dirct(2)*(1-cos(theta))-dirct(3)*sin(theta)   dirct(1)*dirct(3)*(1-cos(theta))+dirct(2)*sin(theta);
    dirct(1)*dirct(2)*(1-cos(theta))+dirct(3)*sin(theta)  dirct(2).^2*(1-cos(theta))+cos(theta)   dirct(2)*dirct(3)*(1-cos(theta))-dirct(1)*sin(theta);
    dirct(1)*dirct(3)*(1-cos(theta))-dirct(2)*cos(theta)  dirct(2)*dirct(3)*(1-cos(theta))+dirct(1)*sin(theta)   dirct(3).^2*(1-cos(theta))+cos(theta);
    ];

origin=repmat(origin(:)',size(P,1),1);
Pr=(P-origin)*T'+origin;