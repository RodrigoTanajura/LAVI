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
