function [M, K] = no_especial(nn, M, K, Ke2)
%Esta função adiciona o valor dos nós especiais às matrizes M e K.
%Carrega os dados do sistema
%Dados_Friswell
%Dados_Eixos
ne = nn - 1;
%Escreve na matriz
J1 = 7e3;
M(1,1) = M(1,1) + J1;
M((ne/2 + 1), (ne/2 + 1)) = M((ne/2 + 1), (ne/2 + 1))*(1 + n^2)...
+ J2 + J3*n^2;
M((ne/2 + 1), (ne/2 + 1) + 1) = M((ne/2 + 1), (ne/2 + 1) + 1)*(-n);
M((ne/2 + 1) + 1, (ne/2 + 1)) = M((ne/2 + 1) + 1, (ne/2 + 1))*(-n);
K((ne/2 + 1), (ne/2 + 1)) = K((ne/2 + 1), (ne/2 + 1)) - Ke2(1,1)...
+ Ke2(1,1)*n^2;
K((ne/2 + 1), (ne/2 + 1) + 1) = K((ne/2 + 1), (ne/2 + 1) + 1)*(-n);
K((ne/2 + 1) + 1, (ne/2 + 1)) = K((ne/2 + 1) + 1, (ne/2 + 1))*(-n);
M(nn, nn) = M(nn, nn) + J4;
end