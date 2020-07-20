%% Complex pressure simultion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Por Sidney Volney Cândido (2020)
% input: nps_low = Limite inferior
        %nps_high = Limite superior
        %fase = pi ou 2pi
        %thetas = vetor de angulos
% output: A = Pressão complexa aleatória
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = complex_pressure(nps_low, nps_high, fase, thetas)
% nps_low = 100;
% nps_high = 140;  
%  equação para ter valores aleatórios entre os limites
nps = nps_low + (nps_high-nps_low).*rand(length(thetas),1);
% converte para pressão     
po = round(db2mag(nps)'); 

% Gerando fases aleatórias
th180 = deg2rad(fase);
fase = th180.*rand(length(thetas),1)';


% Pressão complexa  
A = (po.*exp(1j*fase)); 
