function [EG, f] =  entradas_tcc
% Só uma função que serve para otimizar o espaço da rotina principal
% Gerais
%         EG.df = round((10000-20)/400);
        % Frequência de análise [Hz]
        EG.freq_sim =exp(log(50):0.1:log(10000))'; 
        EG.freq = exp(log(80):.01:log(10000)); 
        
        % Frequência angular [rad/s]
        EG.w = 2*pi*EG.freq_sim; 
        % Temperatura
        EG.To = 20;
        % Pressão estátitca
        EG.Po = 101325; 
        % Umidade relativa
        EG.HR = 78; 
        % Caracteristicas do ar
    [EG.rho0,EG.c0,EG.vis,EG.gam,EG.B2,EG.Cp,EG.kappla]=propair_panneton(EG.To,EG.Po,EG.HR);
        % Número de onda no ar
        EG.k0 = 2*pi.*EG.freq_sim/EG.c0;
        
        % Impedancia caracteristica do ar
        EG.z0 = EG.rho0*EG.c0; 
        % Configura plot
        f = config_plot;
        EG.co = EG.c0;
        
        return

        