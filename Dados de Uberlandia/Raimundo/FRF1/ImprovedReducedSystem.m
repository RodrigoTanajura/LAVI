function [KR,MR,Wirs] = ImprovedReducedSystem(K,M,SlaveDofs)
%Guyan---------------------------------------------------------%
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
%--------------------------------------------------------------%
%Reorganizando as matrizes do sistema original
M = [Mmm, Msm’; Msm, Mss];
K = [Kmm, Ksm’; Ksm, Kss];
%--------------------------------------------------------------%
P = - inv(Kss) * Ksm;
tamanhos = size(P);
W = [ eye(tamanhos(2)); P ];
MR = W’*M*W;
KR = W’*K*W;
%--------------------------------------------------------------%
S = [Kmm*0, Ksm’*0; Ksm*0, inv(Kss)];

Wirs = W + S*M*W*inv(MR)*KR ;
MR = Wirs’*M*Wirs;
KR = Wirs’*K*Wirs;