clc; 
close all; 
clear;
%Carlos Lima Santoro
%Script para comparar diferentes métodos de redução de modelo
%Monta o sistema
%DadosFriswell
%DadosEixos
%Sistema1Eixos
%Aplica as reduções de modelo
[KR(:,:,1), MR(:,:,1), T(:,:,1), Kss] = GuyanReduction(K , M ,SlaveDofs);
[KR(:,:,2), MR(:,:,2), T(:,:,2)] = ImprovedReducedSystem(K , M , SlaveDofs);
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
