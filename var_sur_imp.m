%% Varyng surface impedance on the oblique incidence angle 
% "Assuming that propagation distance of the incident and reflection signal
% inside the material at normal incidence is equal to material thickness
% t... PAN" (Ref)
% Ref: In situ measurements of absorption characteristics using two mics
% and EAnoise (Takahashi et. al)
 
function [zs_field] = var_sur_imp(zc, kc, theta_vector,thick,freq)

    for k = 1:length(theta_vector)
    angulo = theta_vector(k);
    % PAN surface impedance 
    zs_theta(:,k) = zc.*coth(kc.*(thick/cos(angulo)));
    % Reflection Coeficient 
    R(:,k) = (zs_theta(:,k) - 1)./(zs_theta(:,k) + 1);
    % Absorption Coeficient
    abs_theta(:,k) = 1 - abs(R(:,k)).^2;
    end
%% Field incidence (Paris Formula) 

    passo2 = sum(sin(theta_vector).*cos(theta_vector));
    
for g = 1:length(freq)
    % Surface impedance
        passo1 = sum(zs_theta(g,:).*sin(theta_vector).*cos(theta_vector));

        zs_field(g,1) = passo1/passo2;
    

    end

     % Reflection Coeficient 
    R= (zs_field - 1)./(zs_field + 1);
    % Absorption Coeficient
    abs_field = 1 - abs(R).^2;
    


end
