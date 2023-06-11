function [KR,MR,W] = SEREP3(K, M, SlaveDofs)
Dof = length(K(:,1));
SlaveDofs = sort(SlaveDofs);
index = 1 : Dof ;
index(SlaveDofs) = [];
for i = 1 : length(index)
for j = 1 : length(index)
Mmm(i,j) = M(index(i),index(j));
Kmm(i,j) = K(index(i),index(j));
end
end
for i = 1 : length(SlaveDofs)
for j = 1 : length(SlaveDofs)
Mss(i,j) = M(SlaveDofs(i),SlaveDofs(j));
Kss(i,j) = K(SlaveDofs(i),SlaveDofs(j));
end
end
for i = 1 : length(SlaveDofs)
for j = 1 : length(index)
Msm(i,j) = M(SlaveDofs(i),index(j));
Ksm(i,j) = K(SlaveDofs(i),index(j));
end
end
%---------------------------------------------------------------%
%Reorganizando PHIsm matrizes do sistema original
M = [Mmm, Msm'; Msm, Mss];
K = [Kmm, Ksm'; Ksm, Kss];
%---------------------------------------------------------------%
[autovetores,autovalores] = eig(K,M);
tam = length(index);
PHImm = autovetores(1:tam, 1:tam);
PHIsm = autovetores((tam + 1):length(autovetores), 1:tam);
W = [PHImm;PHIsm]*PHImm^(-1);
MR = W'*M*W;
KR = W'*K*W;
