function Z_optmized=dif_imp_optmizer_pp(f,c0,rho0 ,coord,Z_avg,hab, A, step)
% Adapted from brando Z_optmized_pp of insitu_impdance_phd-master
% repository
% Sidney Volney Candido

% Calculates surface impedance by q-term interation process
%Limit the frequency range to speed up the calculation
Range=1:step:length(f);       % step = Every nth frequency points is taken
Range(f(Range)<100)=[];     % Frequencies below 100 Hz are not calculated
Range(f(Range)>10000)=[];   % Frequencies above 10000 Hz are not calculated

%% Geometry

Zf=Z_avg;         % Input impedance
% Variables:
% Here          Name in Finn's paper
% D         :   D'
% Zf(b)     :   Zm (measured impedance at the sensor location)
% Zg        :   Z' or Zk+1 (next surface impedance guess)
% Zg1       :   Zk (guessed surface impedance)
% Zg2       :   Zk-1 (previous guessed surface impedance)
% fZg       :   f(Z') (difference between measured Zm and calucation with Zg)
% fZg1      :   f(Z') (difference between measured Zm and calucation with Zg1)
% fZg2      :   f(Z') (difference between measured Zm and calucation with Zg2)
%%%%%%%%%%%%%%%%%
% hab = transfer function between a and b (b/a)
% A = Complex pressure vector (for sum_pressure2.m)
% coord = struct of coordinates of the experiment geometry (for
% sum_pressure2.m)

% Iniciate
Z_optmized(:,1)=zeros(length(f),1);
hqterm = waitbar(0,'Calculating Z by q-term Model...');
for b=Range
    waitbar((b)/max(Range),hqterm)
    k =2*pi*f(b)/c0; %Wavenumber
    for a=1:40
            if a==1         %first iteration the guessed impedance is the measured impedance
                Zg=Zf(b);   Zg1=1;      fZg=1;
            elseif a==2    %second iteration the guessed impedance is first calculated imp
                            Zg1=Zg;     fZg1=fZg;   Zg=Zg-fZg;
            else            %all other iterations the secant method is used to guess the impedance
                if isfinite(Zg) && Zg~=Zg1 && fZg1~=fZg % To avoid errors
                            Zg2=Zg1;    fZg2=fZg1;
                            Zg1=Zg;     fZg1=fZg;

                    Zg=Zg1-(Zg1-Zg2)/(fZg1-fZg2)*fZg1;
                elseif isinf(Zg) ||  isnan(Zg)
                    break
                end
            end

            if isfinite(Zg) && Zg~=Zg1 && abs(fZg)>0.000001; % To avoid errors         
                
                p1=sum_pressure2(f(b), c0, rho0, Zg, coord.mic1, coord.sph,A);
                
                

                p2=sum_pressure2(f(b), c0, rho0, Zg, coord.mic2, coord.sph, A);
                
                fZg=hab(b)-p2./p1;
            end
            if abs(fZg)<0.000001
                break
            end

    end
    %Zfinn(b,1)=Zg1;
    Z_optmized(b,1)=Zg1;
end
close(hqterm)


