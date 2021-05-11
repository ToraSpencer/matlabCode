function [] = Three_dim_bezier_video()
% Three_dim_bezier_video : 3D Bezier curver construction video with De Casteljau algorithm
%
% Author & support : nicolas.douillet (at) free.fr, 2020.
step_t = 1e-2;
t = 0:step_t:1;
s = 1-t;
% 3D control points coordinates definition
A = [6 0 5]';
B = [0 18 4]'; % 16 2 4
C = [21 25 11]';
D = [33 13 24]';
% De Casteljau algorithm :
% iterative construction
% of M(t) Bezier curve points
E = zeros(3,length(t));
F = zeros(3,length(t));
G = zeros(3,length(t));
H = zeros(3,length(t));
I = zeros(3,length(t));
filename = '3D_bezier_DeCasteljau_construction_algo_bf.gif';
h = figure;
set(h,'Position',[1 1 1920 1080]);
axis tight manual;
for k = 1:length(t)
                
    E(:,k) = s(1,k)*A + t(1,k)*B;
    F(:,k) = s(1,k)*B + t(1,k)*C;
    G(:,k) = s(1,k)*C + t(1,k)*D;
    H(:,k) = s(1,k)*E(:,k) + t(1,k)*F(:,k);
    I(:,k) = s(1,k)*F(:,k) + t(1,k)*G(:,k);
    M(:,k) = s(1,k)*H(:,k) + t(1,k)*I(:,k);
    
    % Display points
    plot3(A(1),A(2),A(3),'co','Linewidth',2), hold on;
    plot3(B(1),B(2),B(3),'co','Linewidth',2), hold on;
    plot3(C(1),C(2),C(3),'co','Linewidth',2), hold on;
    plot3(D(1),D(2),D(3),'co','Linewidth',2), hold on;
    % Display segments
    plot3([A(1),B(1)],[A(2),B(2)],[A(3),B(3)],'Color',[0 1 1],'Linewidth',2), hold on;
    plot3([B(1),C(1)],[B(2),C(2)],[B(3),C(3)],'Color',[0 1 1],'Linewidth',2), hold on;
    plot3([C(1),D(1)],[C(2),D(2)],[C(3),D(3)],'Color',[0 1 1],'Linewidth',2), hold on;
    plot3(E(1,k),E(2,k),E(3,k),'go','Linewidth',2), hold on;
    plot3(F(1,k),F(2,k),F(3,k),'go','Linewidth',2), hold on;
    plot3(G(1,k),G(2,k),G(3,k),'go','Linewidth',2), hold on;
    
    plot3([E(1,k),F(1,k)],[E(2,k),F(2,k)],[E(3,k),F(3,k)], 'Color', 'g', 'Linewidth',2), hold on;
    plot3([F(1,k),G(1,k)],[F(2,k),G(2,k)],[F(3,k),G(3,k)], 'Color', 'g', 'Linewidth',2), hold on;
        
    plot3(H(1,k),H(2,k),H(3,k),'o', 'Color', 'g', 'Linewidth',2), hold on;
    plot3(I(1,k),I(2,k),I(3,k),'o', 'Color', 'g','Linewidth',2), hold on;
    
    plot3([H(1,k),I(1,k)],[H(2,k),I(2,k)],[H(3,k),I(3,k)], 'Color', 'y', 'Linewidth',2), hold on;  
    
    plot3(M(1,:),M(2,:),M(3,:),'r','Linewidth',2);
    plot3(M(1,k),M(2,k),M(3,k),'ro','Linewidth',2);            
    
    set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1], 'FontSize', 14);
    view(-58, 21);
    grid on;
    
    xlabel('X', 'Color', [1 1 1]), ylabel('Y', 'Color', [1 1 1]), zlabel('Z', 'Color', [1 1 1]);
    title('3D Bezier curve built with De Casteljau algorithm', 'Color', [1 1 1], 'FontSize', 16);
    
    drawnow;
    
    frame = getframe(h);    
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the .gif file
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',Inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    
    clf(1);        
end
end