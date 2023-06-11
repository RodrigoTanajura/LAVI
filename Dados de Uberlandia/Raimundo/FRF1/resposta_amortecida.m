function [omega,X0] = resposta_amortecida(M,K,F,intervalo,d_omega)
%Funcao para calcular a resposta em frequencia de um sistema amortecido
% Problema: (K + w*C*j - w^2M)*X0 = F
aux = 1;
for i = intervalo(1) : d_omega : intervalo(2)
omega(aux) = i;
X0(aux,:) = (K - (omega(aux)^2)*M)\F;
aux = aux+1;
end
end