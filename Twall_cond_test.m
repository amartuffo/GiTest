function [Twall, Q] = Twall_cond_test( H, rhoCpS, ks, q_IPS, Twall, icetype, Twat)

% Twall è il VETTORE che contiene la distribuzione di temperatura nella
% parete 

%DA CONTROLLARE (soprattutto Twater) secondo me la Twater non è da mettere
%se non al primo step. Successivamente infatti viene presa come Twater a
%h=0 la Twall(end) del passo precedente. Inoltre non bisogna modificare le
%altre funzioni perché esse si basano già sulla Twall(end).

global Beta LWC r hAir hWater Pinf Tinf Vinf cpAir cpWater dt 

%% Determinazione coefficiente convettivo tra substrato e ghiaccio
% Bisogna definire il valore di h_conv in base al "tipo di ghiaccio" che si
% forma sopra al substrato.

% 0-AS, 1-AWS, 2-AIS, 3-AWIS, 4-AIWS, 5-AWIWS
if icetype == 0 
    h_conv = hAir;            %hAir = 1000 [W/(K*m^2)]
    
elseif icetype == 1
    h_conv = 30000;          %hWater = 30000 [W/(K*m^2)] Moving water film
    
elseif icetype == 2 || icetype == 3
    h_conv = 1000;            %hIceS = ??? [W/(K*m^2)]
    
elseif icetype == 4 || icetype == 5
    h_conv = 1000;            %hWater_static = ??? [W/(K*m^2)] Static water film
    
end    


%% Discretizzazione substrato

N = 50; %(Può essere modificato in base al numero di cv) 
dz = H/N;
z = [0 dz/2:dz:H-dz/2 H];


%% Definizione temperature iniziali

% Se non viene specificata la T dell'acqua all'interfaccia allora assumiamo
% che sia uguale a Twall all'istante n 
if nargin == 6
    Twat = Twall(end);
end

T0 = Twall;            % Chiamiamo T0 la Twall all'istante n
T1 = Twall;               % Per iniziare, assegnamo Twall(n+1)=Twall(n) 


%% Modello Conduzione

err = 1;
eps = 1e-6;

while err > eps
    
    Tstim = T1; % Stimatore di Twall(n+1) per il calcolo dell'errore
    
    % Nodi interni
    for i = 2 : N+1
        
        if i == 2 
            a = ((dz*rhoCpS)/dt)+3*ks/dz;
            b = 2*ks/dz;    
            c = ks/dz;
            d = ((dz*rhoCpS)/dt)*T0(i) + q_IPS/N;  %divido per il numero di cv
        elseif i == N+1
            a = ((dz*rhoCpS)/dt)+3*ks/dz;
            b = ks/dz;    
            c = 2*ks/dz;
            d = ((dz*rhoCpS)/dt)*T0(i) + q_IPS/N;  %divido per il numero di cv
        else
            a = ((dz*rhoCpS)/dt)+2*ks/dz;
            b = ks/dz;    
            c = ks/dz;
            d = ((dz*rhoCpS)/dt)*T0(i) + q_IPS/N;  %divido per il numero di cv
        end
        
        T1(i) = (b*T1(i-1)+c*T1(i+1)+d)/a;
    
     end
    
     % Primo Nodo (Hp: parete adiabatica)
     T1(1) = T1(2);  
   
     % Ultimo Nodo (Considero convezione con lo strato di fluido superiore)
     if icetype == 0
                  
         %chiEvap = 11*6.3e4/Pinf;
         %e0      = 27.03;
         q_evap = EvapHeatFlux(T1(end),Tinf);
         chiEvap = q_evap/(T1(end)-Tinf);
         e0 = 1;
         
         q_aero  = 0.5 * r * hAir * Vinf^2 / cpAir;
         q_kin   = 0.5 * LWC * Beta * Vinf^3;
         
         a = 2*ks/dz + h_conv + Beta * LWC * Vinf * cpWater + chiEvap * e0;
         b = 2*ks/dz;
         c = (h_conv + chiEvap * e0 + Beta * LWC * Vinf * cpWater)*Tinf + q_aero + q_kin;
         
         T1(N+2) = (b*T1(N+1)+c)/a;
     else
         
         a = 2*ks/dz + h_conv;
         b = 2*ks/dz;
         c = h_conv*Twat;  
         
         T1(N+2) = (b*T1(N+1)+c)/a;
         
     end

     err = max(abs(Tstim-T1));

end

if icetype == 0
         q_aero  = 0.5 * r * hAir * Vinf^2 / cpAir
         q_kin   = 0.5 * LWC * Beta * Vinf^3
         q_drop  = Beta * LWC * Vinf * cpWater * (T1(end) - Tinf)
         q_evap  = EvapHeatFlux(T1(end),Tinf)
         q_conv_AS = h_conv*(T1(end) - Tinf)
         
else
    q_conv_AWS = h_conv*(T1(end) - Twat)
end

Twall = T1; 
Q = - 2* ks * (T1(end) - T1(end-1))/dz;



