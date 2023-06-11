clc; close all; clear;
%Carlos Lima Santoro
%Script para comparar diferentes métodos de redução de modelo
%Monta o sistema
%DadosFriswell
%DadosEixos
%Sistema1Eixos
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
[M, K] = noespecial(nn, M, K, Ke2);

%Definição dos graus de liberdade escravos
SlaveDofs = 2:(nn-1);
%Definição dos vetores de forçamento
F = [1; zeros(nn - 1, 1)];
FR = [1; zeros(nn - length(SlaveDofs) - 1, 1)];
function [M, K] = no_especial(nn, M, K, Ke2)
%Esta função adiciona o valor dos nós especiais às matrizes M e K.
%Carrega os dados do sistema
%Dados_Friswell
%Dados_Eixos

ne = nn - 1;
%Escreve na matriz
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

%Aplica as reduções de modelo
[KR(:,:,1), MR(:,:,1), T(:,:,1), Kss] = GuyanReduction(K , M ,SlaveDofs);
[KR(:,:,2), MR(:,:,2), T(:,:,2)] = ImprovedReducedSystem(K , M ,SlaveDofs);
[KR(:,:,3), MR(:,:,3), T(:,:,3)] = SEREP(K , M , SlaveDofs);
%Calcula a resposta do sistema completo
intervalo = [0, 1000];
d_omega = 0.1;
[omega,x0] = calcularesposta(M,K,F,intervalo,d_omega);
%Calcula a resposta dos sistemas reduzidos
for i = 1:(length(KR))
    [omega,x0R(:,:,i)] = calcularesposta(MR(:,:,i),KR(:,:,i),FR,intervalo,d_omega);
end
%Coloca tudo em um vetor para ser plotado no mesmo gráfico
X(:,1) = x0(:,1);
X(:, 2:(length(KR) + 1)) = x0R(:, 1, 1:length(KR));

%Plota os resultados
plot(omega,20*log10(abs(X)))
nomes = ["Completo", "Guyan", "IRS", "SEREP"];
legend(nomes)
xlabel('Frequência (rad/s)')
ylabel('Amplitude (dB)')
title('FRF Rotor')

function [KR,MR,W,Kss]=GuyanReduction(K , M , SlaveDofs)

% Redução de Guyan

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
%Reorganizando as matrizes do sistema original
M = [Mmm, Msm'; Msm, Mss];
K = [Kmm, Ksm'; Ksm, Kss];
%---------------------------------------------------------------%
P = - inv(Kss) * Ksm;
tamanhos = size(P);
W = [ eye(tamanhos(2)); P ];
MR = W'*M*W;
KR = W'*K*W;
end

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
M = [Mmm, Msm'; Msm, Mss];
K = [Kmm, Ksm'; Ksm, Kss];
%--------------------------------------------------------------%
P = - inv(Kss) * Ksm;
tamanhos = size(P);
W = [ eye(tamanhos(2)); P ];
MR = W'*M*W;
KR = W'*K*W;
%--------------------------------------------------------------%
S = [Kmm*0, Ksm'*0; Ksm*0, inv(Kss)];
Wirs = W + S*M*W*inv(MR)*KR ;
MR = Wirs'*M*Wirs;
KR = Wirs'*K*Wirs;
end

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
end

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

clc; close all; clear;
%Carlos Lima Santoro
%Script para comparar diferentes métodos de redução de modelo
%Monta o sistema
DadosFriswell
DadosEixos
Sistema1Eixos
%Aplica as reduções de modelo
[KR(:,:,1), MR(:,:,1), T(:,:,1), Kss] = GuyanReduction(K , M,...
SlaveDofs);
[KR(:,:,2), MR(:,:,2), T(:,:,2)] = ImprovedReducedSystem(K , M,...
SlaveDofs);
[KR(:,:,3), MR(:,:,3), T(:,:,3)] = SEREP(K , M , SlaveDofs);
%Calcula a resposta do sistema completo
intervalo = [0, 1000];
d_omega = 0.1;
[omega,x0] = calcularesposta(M,K,F,intervalo,d_omega);
%Calcula a resposta dos sistemas reduzidos
for i = 1:(length(KR))
    [omega,x0R(:,:,i)] = calcularesposta(MR(:,:,i),KR(:,:,i),...
    FR,intervalo,d_omega);
end
%Coloca tudo em um vetor para ser plotado no mesmo gráfico
X(:,1) = x0(:,1);
X(:, 2:(length(KR) + 1)) = x0R(:, 1, 1:length(KR));
%Calcula os Dr e orgniza os graus de liberdade em seus respectivos...

lugares
tam = size(x0R);
for i = 1:tam(3)
    Dr(:,:,i) = T(:,:,i)*x0R(:,:,i)';
    aux = Dr(2,:,i);
    tam = size(Dr);
    for j = 3:(tam(1) - 1)
        Dr(j,:,i) = Dr(j + 1, :, i);
    end
    Dr(tam(1),:,i) = aux;
end
%Aplicação da euqação 11 do artigo do Sellgren
tam = size(Dr);
for j = 1:tam(3)
    for i = 1:length(Dr(:,:,1))
        gamma(i,j) = 1 - abs(Dr(:,i,j)'*x0(i,:)'/(norm(Dr(:,i,j))*norm(x0(i,:))));
     end
end
plot1 = plot(omega, gamma)
set(plot1(1), 'Color', [1 0 0])
set(plot1(2), 'Color', [0 1 0])
set(plot1(3), 'Color', [0 0 1])
nomes = ["Guyan", "IRS", "SEREP"];
legend(nomes)
xlabel('FREQUÊNCIA')
ylabel('ERRO RELATIVO')
end


