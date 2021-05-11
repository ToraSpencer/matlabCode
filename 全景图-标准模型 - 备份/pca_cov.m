function [C,M] = pca_cov(vertex)
C = mean(vertex);
P = vertex-C;
M_cov = [sum(P(:,1).*P(:,1))/(length(P(:,1))-1) sum(P(:,1).*P(:,2))/(length(P(:,1))-1) sum(P(:,1).*P(:,3))/(length(P(:,1))-1);...
         sum(P(:,2).*P(:,1))/(length(P(:,1))-1) sum(P(:,2).*P(:,2))/(length(P(:,1))-1) sum(P(:,2).*P(:,3))/(length(P(:,1))-1);...
         sum(P(:,3).*P(:,1))/(length(P(:,1))-1) sum(P(:,3).*P(:,2))/(length(P(:,1))-1) sum(P(:,3).*P(:,3))/(length(P(:,1))-1)];
[M,~] = eig(M_cov);
end