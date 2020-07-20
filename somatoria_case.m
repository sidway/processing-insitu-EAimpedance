function [z_avg, alpha, hab] = somatoria_case(input, coord, z_input)
%% 3. Presão complexa
randp.A = complex_pressure(20, 120, 2*pi, coord.thetas);
%% 4. Simulação da pressão
% Mic 1
    pressure.pa = sum_pressure2(input.freq_sim, input.c0, input.rho0,...
        z_input*(input.rho0*input.c0), coord.mic1, coord.sph, randp.A);
% Mic 2 (PP e PU)
    pressure.pb = sum_pressure2(input.freq_sim, input.c0, input.rho0,...
        z_input*(input.rho0*input.c0),  coord.mic2, coord.sph, randp.A);
%% 5. Estimativas da impedância de superfície
   
hab = pressure.pb./ pressure.pa; % FT entre os microfones 
    l = coord.mic1(3) - coord.mic2(3);
    d = coord.mic2(3);
[z_avg] = ra_pp_estimation(hab,input.k0, l, d,input.rho0,input.c0); % Estimativa z_avg
[~, alpha] = reflection_and_absorption_coefficient(z_avg,input.z0); % Absorção

