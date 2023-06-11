clc; close all; clear;
%Carlos Lima Santoro
%Script para comparar diferentes métodos de redução de modelo
%Monta o sistema
%DadosFriswell
%DadosEixos
%Sistema1Eixos
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
%Calcula os Dr e orgniza os graus de liberdade em seus respectivos...lugares
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