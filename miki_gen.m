%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% === Surface impedance calculation with MIKI coeficients === %%%
% by Sidney Volney Cândido in April 2020
%
% Ref: Y. MIKI - Modifications of Delany-Bazley Models (1990)
%
% Inputs:
%           resf = Flow resistivity
%           f = Frequency vector;
%           ro = Air density
%           co = Speed of sound in air
%           d1 = Sample thickness
% Outputs:  
%           kc = Propagation Constant (refered as y in the paper)
%           zc = Characteristic impedance (refered as z0 in the paper)
%           zs = surface impedance for normal incidence
%           vp = Absorption coeficient
%% =========================================================  %%
function [kc, zc, zs, vp] = miki_gen(resf, f, ro, co, d1)
         
% Demo (coment function)
%     resf=25000; %Resistividade ao fluxo [kg*s^-1*m^-3]
%     f=100:5:10000; %Vetor Frequencia [Hz]
%     ro=1.2; %Densidade do ar [kg/m³]
%     co=343; %Velocidade do ar [m/s]
%      d1=[25*10^-3]; 
    % Angular Frequency
    w=2*pi.*f;
    % Wavenumber in air
    k0=(2*pi/co).*f;
    % Ration Frequency flow resitivity
    V=(f/resf);
    
    % Impedancia caracteristica normalizada
            stp1 = 1+0.07*(V).^-.632;
            stp2 = -.107*(V).^-.632;
            zc = (stp1 + 1j*stp2);
            
    % Costante de propagação
            stg1 = (w/co).*(0.16*(V).^-.618);%;
            stg2= (w/co).*(1+ 0.109*(V).^-.618);%;
            kc = (stg1 + 1j*stg2);%.*k0;
            
            zs=zc.*coth(kc.*d1); %Impedancia de superficie
            vp=1-(abs((zs-1)./(zs+1))).^2; %Coeficiente de absorção

return