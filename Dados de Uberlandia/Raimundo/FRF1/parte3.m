N2 = 56;
N3 = 111; %numero de dentes
n = N2/N3; %razão de engrenamento
J1 = 7e3;
J2 = 4e3;
J3 = 13e3;
J4 = 5e3; %Inércias [kgm^2]
k1 = 1.6e9;
k2 = 0.29e9; %Rigidez dos eixos [Nm/rad]
%Propriedades gerais do aço
ro = 7860; %Densidade do aço em kg/m^3
G = 80e9; %Módulo de cisalhamento em Pa
l1 = 1; %Comprimento do eixo 1 em metros
l2 = 1; %Comprimento do eixo 2 em metros
%Diâmetros em função de k, l e G
d1 = (32*k1*l1/(G*pi))^(1/4); %Diâmetro do eixo 1 em metros
d2 = (32*k2*l2/(G*pi))^(1/4); %Diâmetro do eixo 2 em metros
m_eixo1 = ro*l1*pi*d1^2/4; %Massa do eixo 1 em kg
m_eixo2 = ro*l2*pi*d2^2/4; %Massa do eixo 2 em kg
Je1 = m_eixo1*(d1/2)^2/2; %Inércias de giro dos eixos em kg.m^2
Je2 = m_eixo2*(d2/2)^2/2;
disp(Je1/J1);
%Esta rotina define o sistema com base nos dados de Friswell e nos dados
%dos eixos, discretizando o sistema em diversos pedaços.
nepe = 5; %Número de elementos por eixo
ne = 2*nepe; %Número de elementos
nn = ne + 1; %Número de nós
le1 = l1/nepe; %Comprimento do elemento do eixo 1
le2 = l2/nepe; %Comprimento do elemento do eixo 2
A1 = pi*d1^4/32; %Momento polar de área do eixo 1
A2 = pi*d2^4/32; %Momento polar de área do eixo 2
Ke1 = (G*A1/le1)*[1,-1;-1,1]; %Matriz de rigidez dos elementos do eixo 1
Ke2 = (G*A2/le2)*[1,-1;-1,1]; %Matriz de rigidez dos elementos do eixo 2
Me1 = (ro*A1*le1/6)*[2,1;1,2]; %Matriz de massa dos elementos do eixo 1
Me2 = (ro*A2*le2/6)*[2,1;1,2]; %Matriz de massa dos elementos do eixo 2
nos = zeros(ne,2); %Matriz de nós (cada posição da matriz é relativa...a um elemento)
for i = 1:ne
nos(i, 1:2) = [i, i+1];
end
%Montando as matrizes de massa e rigidez
M = zeros(nn,nn);
K = M;
for i = 1:(ne/2)
M(nos(i,1):nos(i,2), nos(i,1):nos(i,2)) = M(nos(i,1):nos(i,2),...
nos(i,1):nos(i,2)) + Me1;
K(nos(i,1):nos(i,2), nos(i,1):nos(i,2)) = K(nos(i,1):nos(i,2),...
nos(i,1):nos(i,2)) + Ke1;
end
for i = (ne/2 + 1):ne
M(nos(i,1):nos(i,2), nos(i,1):nos(i,2)) = M(nos(i,1):nos(i,2),...
nos(i,1):nos(i,2)) + Me2;
K(nos(i,1):nos(i,2), nos(i,1):nos(i,2)) = K(nos(i,1):nos(i,2),...
nos(i,1):nos(i,2)) + Ke2;
end
[M, K] = no_especial(nn, M, K, Ke2);
%Definição dos graus de liberdade escravos
SlaveDofs = 2:(nn-1);
%Definição dos vetores de forçamento
F = [1; zeros(nn - 1, 1)];
FR = [1; zeros(nn - length(SlaveDofs) - 1, 1)];