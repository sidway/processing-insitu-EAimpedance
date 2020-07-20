
%% Parte da function
function pressure = sum_pressure2(frequencia, c0, rho0, zs, coord_mic, com_sph, A)
    porous.Z = zs; %%% impedância de entrada
     %%% Vetor de frequência
    vel_air = c0; %%% Velocidade do som
    den_air = rho0; %%% Desindade do ar
    
    % Coordenadas do microfones no plano xy
    cx = coord_mic(1);
    cy = coord_mic(2);
    
     % Posicação z da fonte
    hs = com_sph(3,:);
     % Posição z do receptor
    za = coord_mic(3);
    % Distancia horizontal da fonte
    r = sqrt((com_sph(1,:) - cx).^2 + (com_sph(2,:) - cy).^2);

    % Distancia da fonte real e imagem até o sensor
     R1=sqrt((r.^2)+((hs-za).^2));
     R2=sqrt((r.^2)+((hs+za).^2));
     

% pressure = zeros(1, length(frequencia));

    beta = (den_air*vel_air) ./ porous.Z;
    k0 = 2 * pi * frequencia/vel_air;
    
    % Ponto A (distancia l+d)
    % Integrand
    int_pa=@(q)(sum(A.*((exp((-q.*k0).*beta)).*...
        exp(-1i*k0.*(sqrt((r.^2)+(hs+za-1i*q).^2))))./...
        (sqrt((r.^2)+(hs+za-1i*q).^2)),2));
    % Integral
    Iqa=-(2*k0.*beta).*integral(int_pa,0,20,'ArrayValued',true);
    pressure=sum(A.*(exp(-1i*k0*R1)/R1),2)+...
        sum(A.*(exp(-1i*k0*R2)/R2),2)+sum(Iqa,2);



% end
