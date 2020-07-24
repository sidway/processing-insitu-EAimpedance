function resist_fit=bazley_alpha_fit2(xdata,ydata,dimm,x0)
%%% This function uses the least square method to fit experimental data
%%% with Bazley impedance model. It finds the correct flow resistivity

%%% Input parameters:
%%% f: Frequency vector of the measurement
%%% alpha: absorption vector (same size of f)
%%% c0: sound speed
%%% rho0: fluid density
%%% dimm: sample's thickness in [mm]
%%% x0: initial guess of flow resistivity

%% Function to fit

d1=dimm/1000;     % Thickness in [mm]

zs_field = @(resist,xdata) material_reference2(xdata, d1, resist);
R=@(resist,xdata)((zs_field(resist,xdata)-1)./(zs_field(resist,xdata)+1));
alpha=@(resist,xdata)(1-abs((R(resist,xdata)).^2));
%% Curve fitting
resist_fit=lsqcurvefit(alpha,x0,xdata,ydata);