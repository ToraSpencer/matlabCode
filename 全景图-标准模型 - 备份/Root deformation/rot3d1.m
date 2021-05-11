function Pr=rot3d1(P,origin,dirct,theta)
% �������P���ţ���origin�㣬����Ϊdirct��ֱ�ߣ���תtheta��
% P����Ҫ��ת�����꼯�ϣ�n��3����
% origin��ת��ͨ���ĵ㣬1��3����
% direct��ת�᷽��������1��3����
% theta����ת�Ƕȣ���λ����


dirct=dirct(:)/norm(dirct);
T = [dirct(1).^2*(1-cos(theta))+cos(theta)  dirct(1)*dirct(2)*(1-cos(theta))-dirct(3)*sin(theta)   dirct(1)*dirct(3)*(1-cos(theta))+dirct(2)*sin(theta);
    dirct(1)*dirct(2)*(1-cos(theta))+dirct(3)*sin(theta)  dirct(2).^2*(1-cos(theta))+cos(theta)   dirct(2)*dirct(3)*(1-cos(theta))-dirct(1)*sin(theta);
    dirct(1)*dirct(3)*(1-cos(theta))-dirct(2)*cos(theta)  dirct(2)*dirct(3)*(1-cos(theta))+dirct(1)*sin(theta)   dirct(3).^2*(1-cos(theta))+cos(theta);
    ];

origin=repmat(origin(:)',size(P,1),1);
Pr=(P-origin)*T'+origin;