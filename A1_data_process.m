                %% Processamento de imped�ncia in-situ
                            %%%%%%%%%%%%%%%%%%%%%
                            % Sidney V. C�ndido %
                            %        2020       %
                            %%%%%%%%%%%%%%%%%%%%%
clear; clc;
% Esta � uma rotina mestre. Antes de utilizar, minimize todas as se��es do
% c�digo.
% O objetivo desta rotina � validar o processamento da medi��o de
% imped�ncia in-situ de um material conhecido.
%
% Refer�ncias:
% [1] An exact Laplace transform formulation for a point source above 
% a ground surface. (Di, Xiao and Gilbert, Kenneth E)
% [2] In situ measurements of surface impedance and absorption coefficients
% of porous materials using two microphones and ambient noise (Takahashi,
%  Yasuo and Otsuru, Toru and Tomiku, Reiji)
% [3] A review of the in situ impedance and sound absorption measurement
% techniques (Brand�o, Eric and Lenzi, Arcanjo and Paul, Stephan)

% ATEN��O!!! ESTA ROTINA NECESSITA DO ITATOOLBOX e todas as fun��es do
% pacote
%% 1. Entradas
% Configura��es de plot e vetores de frequencia, numero de onda
[input, f] =  entradas_tcc;
f.fmax = 10000;
f.fmin = 100;

% PrhoPRIEDADES DO AR
    To=18; %[Celsius]
    Po=101100; % [Pa]
    HR=71;  %[por cento]

% Calculo da velocidade do Som e densidade do ar
    [rho,co,~,~,~,~,~]=propair_panneton(To,Po,HR); clear To Po HR
% Experimento
    l = 13/1000;            % Distancia entre micrhofones, intriseco do aparato
    d = 7/1000;             % distancia do mic 1 at� a amostra  % distancia do 
    coord.r = 50;           % Distancia da fonte at� o centro da amostra
    coord.mic1 = [0 0 l+d]; % Coordenada do mic 1
    coord.mic2 = [0 0 d];   % Coordenada do mic 2
    coord.thetas = deg2rad(0:0.25:75);  % vetor de angulos (p/ sum_pressure.m)
    h = length(coord.thetas)+1;        
    coord.sph(2,:) = sin(coord.thetas)*coord.r; % Coordenadas y da fonte
    coord.sph(3,:) = cos(coord.thetas)*coord.r; % Coordenadas z da fonte
% Processamento
    Fs = 51200;             % amostragem na NI
    nfft = 8192;            % Numero de pontos da FFT
    avg = 15;               % Numero de m�dias
    jan = nfft/(avg*0.5);   % Tamanho da janela (nfft/n_med*overlap)
    ks = 15;                % Numero de janelas (parametro do smooth)
% Vetor de amplitude complexa para cada posi��o de fonte
    A = complex_pressure(20, 120, 2*pi, coord.thetas); 
%% 2. Material

thickness = 0.020;  % Espessura do material
resistivity = 12500;% Resistividade ao fluxo

coord_thetas = deg2rad(0:78); % Vetor de angulos (p/ estimativa de 
% incidencia difusa)

% Imped�ncia incidencia difusa (escolher 1)
[z_input] =  material_reference2(input.freq_sim, thickness, resistivity);

% Calculo da absor��o 
[~, alpha_input] = reflection_and_absorption_coefficient(z_input,1);

% Simula��o do experimento (adapta��o de di e gilbert)
% PS: Sem otimiza��o que tira a colora��o 
[zs_field, alpha_field, FT] = somatoria_case(input, coord, z_input);

% Absor��o em ter�o de oitava
[freq_t,alpha_input_third] = narrow_to_one_third_octave(input.freq_sim,alpha_input);
%% 3. Carregando medi��o
cd('C:\TCC\medi��es\2020\med_03_19')
load('med_pp_5.mat')
cd('C:\TCC\Organ_codes\main')
material = 'melamina_preta';
%% 4. Processando medi��es
% tic 
% FT da medi��o normal
[hab_n ,freq] = tfestimate(meas.pressao(:,2), meas.pressao(:,1),jan,[],...
    nfft,Fs);
% FT da medi��o com mics invertidos
[hab_inv ,~] = tfestimate(meas.pressao_inv(:,1), meas.pressao_inv(:,2),...
    jan,[],nfft,Fs);

figure('name', 'Transfer Function');
set(gcf, 'Position', get(0, 'Screensize'));
plot_transferfunction(freq, smooth(hab_n,ks), f, '-*', [0.2 0.0 0.8]); hold on
plot_transferfunction(freq, smooth(hab_inv,ks), f, '--', [0.4 0.2 0.2]); hold on
plot_transferfunction(freq, smooth((hab_n+hab_inv)/2,ks), f, '-', [0.9 0.8 0.2]); hold on
plot_transferfunction(freq, meas.Hm_smooth, f, '--*', [0.4 0.8 0.2]); hold on
%% 5. FILTRO de ter�o de oitava utilizando ITATOOLBOX
% Aqui ser� obtido os valores m�ximos de cada banda de ter�o de oitava
% da FT m�dia para acelerar o processamento da otimiza��o
% Est� em revis�o

hab = itaAudio; % inicia objeto
hab.samplingRate = Fs; % Altera frequencia de amostragem
hab.freq = smooth((hab_n+hab_inv)/2,ks); % Insere a FT no objeto

% Aplica o m�todo
filtro = ita_filter_fractional_octavebands(hab,'bandsperoctave', 3,'freqRange',[100 10000]);
 
oi = filtro.freq(:,:); % Tira do objeto

% Vetor de ter�os de oitava
freq_third = [100 125 160 200 250 ...
                    315 400 500 630 800 1000 1250 1600 2000 2500 ...
                    3150 4000 5000 6300 8000 10000];
% Inicia o vetor de FT por ter�o de oitava
hab_third = zeros(1,length(freq_third));
% Loop que pega os valores maximos
for ms = 1:length(freq_third)
    hab_third(ms) = max(oi(:,ms));
end
%% 6. Dedu��o da imped�ncia de superficie ( Onda plana, campo reverberante )
% Banda Estreita (p/ compara��o)
    % Numero de onda
    k0 = 2*pi*freq_third/co;
    % Dedu��o da imped�ncia de superficie (ondas planas, ruido ambiente)
    [z_avg_third] = ra_pp_estimation(hab_third,k0, l, d,rho,co);
    % Absor��o
    [~, alphab_third] = reflection_and_absorption_coefficient(z_avg_third,...
    rho*co);

% Banda Estreita (p/ compara��o)
    % Numero de onda o vetor de frequencia processado no tfestimate
    k0 = 2*pi*freq/co;

    % M�dia da FT (corre��o de fase
    hm = (hab_n+hab_inv)/2;

    % Dedu��o da imped�ncia de superficie (ondas planas, ruido ambiente)
    [z_avg] = ra_pp_estimation(smooth(hm,ks),k0, l, d,rho,co);
    % Absor��o
    [~, alphab] = reflection_and_absorption_coefficient(z_avg,rho*co);
%% 7. Otimiza��o ( Em ter�o de oitava )
% Imped�ncia de superf�cie
Zq=dif_imp_optmizer_pp(freq_third,co,rho, coord,z_avg_third,hab_third, A,1);
% absor��o 
alpha_q = 1 - (abs((Zq - rho*co)./(Zq + rho*co))).^2;

% Ajustando o plot
Zq(Zq==0)=nan;
Zq(imag(Zq)==0)=1j*nan;
alpha_q(alpha_q==0) = nan;
%% 8. Fit para encontrar resistividade %EM PROGRESSO

resist_fit=difuse_alpha_fit(freq_third(5:14)',alphab_third(5:14)',thickness*1000,8000);
[z_fit] =  material_reference2(freq, thickness, resist_fit);
[~, alpha_fit] = reflection_and_absorption_coefficient(z_fit,1);

figure()
plot_absorption(freq, alpha_fit, f,'--',[1 0.2 0.2]); hold on
plot_absorption(freq_third, alphab_third, f,'--o',[0.7 0.7 0.7]);
plot_absorption(freq_t, alpha_input_third,f,'-',[0.0 0.0 0.0]);
arg = 0;
%% 9. Plots
if arg == 1
%% Impedancia PLOT
cd('C:\TCC\Matriz\Img\all_materiais')
figure('Name', 'Impedancia de superficie')
set(gcf, 'Position', get(0, 'Screensize')); 
plot_impedance(freq, z_avg/(rho*co), f, '--',[0.2 0.2 0.2]);
plot_impedance(freq_third, z_avg_third/(rho*co), f, '--o',[0.7 0.7 0.7]);
plot_impedance(freq_third, Zq/(rho*co), f, '--*',[0.45 0.45 0.45]);
plot_impedance(input.freq_sim, zs_field/(rho*co), f , '-',[0.0 0.0 0.0])


legend('Experimental','Experimental em ter�o de oitava', 'Otimizado',  'Modelo de entrada','location', 'best');%'Otimizado do experimental',
title(['Imped�ncia de superf�cie. n m�dias: ' num2str(avg)...
        '; k smooth = ' num2str(ks)])
    
    saveas(gcf,[material '_impedancia.png'])
    close all
%% Absor��o PLOT
figure('name', 'Absorption Plot');
set(gcf, 'Position', get(0, 'Screensize'));

% Experimental
plot_absorption(freq, alphab, f,'--',[0.2 0.2 0.2]); hold on
% Experimental em ter�o de oitava
plot_absorption(freq_third, alphab_third, f,'--o',[0.7 0.7 0.7]);
% Otimizado 
plot_absorption(freq_third, alpha_q, f,'--*',[0.45 0.45 0.45]);
% Modelo de entrada (referencia)
plot_absorption(freq_t, alpha_input_third,f,'-',[0.0 0.0 0.0]);

legend('Experimental','Experimental ter�o de oitava', 'Otimizado', 'Modelo de entrada','location', 'best');
title(['Coeficiente de absor��o. n m�dias: ' num2str(avg)...
        '; k smooth = ' num2str(ks)])
    
%     saveas(gcf,[material '_absorption.png'])
%     close all
%% Plot Deviation
figure('name', 'Deviation');
set(gcf, 'Position', get(0, 'Screensize'));
plot_absorption(freq_third, abs(alpha_input_third(9:29)-alphab_third),f,'-',[0.0 0.0 0.0]); hold on
plot_absorption(freq_third, abs(alpha_input_third(9:29)-alpha_q'),f,'--o',[0.5 0.5 0.5]); hold on

legend('Experimental','Otimizado','location', 'best');
title(['Desvio abs do coef. de absor��o. n m�dias: ' num2str(avg)...
        '; k smooth = ' num2str(ks)])
    
    saveas(gcf,[material '_deviation.png'])
    close all
cd('C:\TCC\Organ_codes\main')
end