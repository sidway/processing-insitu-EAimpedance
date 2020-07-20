%% Seção Dois Impedância de superficie dependente do campo
function zs_field =  material_reference2(freq, thickness, resitivity)
%% Comparando Miki field com os efeitos anisotropicos
% Entradas Gerais
    [EG, ~] =  entradas_tcc;
    ro = 1.2; co = 343;
    % Material
    EG.t = thickness; EG.resf = resitivity;
 
% Miki Model (incidencia normal)
[kc, zc, zs_normal, ~] = miki_gen(EG.resf, freq, ro, co, EG.t);
theta_vector = deg2rad(0:78);
%% MIKI AVG
zs_field = var_sur_imp(zc, kc, theta_vector,thickness,freq);


end