%% Complex pressure simultion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Por Sidney Volney C�ndido (2020)
% input: nps_low = Limite inferior
        %nps_high = Limite superior
        %fase = pi ou 2pi
        %thetas = vetor de angulos
% output: A = Press�o complexa aleat�ria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = complex_pressure(nps_low, nps_high, fase, thetas)
% nps_low = 100;
% nps_high = 140;  
%  equa��o para ter valores aleat�rios entre os limites
nps = nps_low + (nps_high-nps_low).*rand(length(thetas),1);
% converte para press�o     
po = round(db2mag(nps)'); 

% Gerando fases aleat�rias
th180 = deg2rad(fase);
fase = th180.*rand(length(thetas),1)';


% Press�o complexa  
A = (po.*exp(1j*fase)); 
