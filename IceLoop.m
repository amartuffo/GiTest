clear all 
close all
clc

addpath('\Users\aless\Desktop\IceLoop\IceLoop Definitivo');

% WallTemp  = [274.15, 262.15];
% FarTemp = [272.15, 262.15];
WallTemp  = 1 +273.15;
FarTemp = 1 +273.15;
BoundaryCondition = 1;
fM = 25;
l = 2.5;
pathSave = '/Users/mariachiaragallia/Dropbox/Documents/Università/Progress/Plot_Ice/NewIcingModel/';
save     = false;

PrintValuesIt   = false;
PrintValuesEnd  = true; 
PrintErrorIt    = true; 
PrintErrorEnd   = true; 

PlotTempHeights = false; 
PlotTotal       = true;

IPS             = true;
t_IPS           = 5;
IPS_Heat_Flux   = 20*10^3;      % [W/m^2]      (valori plausibili da 8 a 30)

%Myers Beta = 0.5, LWC = 0.001, Vinf = 90

% fS = figure('Position', get(0, 'Screensize'));

% Two loops, one on the wall temperature to be assigned, one on the
% farfield temperature 

 
for k = 1:length(FarTemp)
    for i = 1:length(WallTemp)
    %% Define simulation parameter
    global Beta LWC Pr r hAir hWater Twall rhoIceRime rhoIceGlaze rhoWater rhoS Tinf Vinf Pinf Tfreezing Lambda_e Lambda_f Lambda_s cpAir cpWater cpS kWater kIce t0 tEnd dt q_IPS;
    global h0 B0 c0 h B c Twater Tdrop Tice TwaterMelt IceRate H;
    global m_imp m_evap m_sub m_water m_ice m_melt m_waterMelt
    global BC 
    global T_interfaccia %aggiunta da me 
    
    BC          = BoundaryCondition;
    Beta        = 0.5;
    LWC         = 0.001;           % [Kg/m^3]
    Pr          = 0.71;
    r           = Pr^(1/3);         
    hAir        = 1000;             % [W/(K*m^2)]
    hWater      = 30000;            % [W/(K*m^2)]
    Twall       = WallTemp(i);      % [K]
    rhoIceRime  = 880;              % [Kg/m^3]
    rhoIceGlaze = 917;              % [Kg/m^3]
    rhoWater    = 1000;             % [Kg/m^3]
    rhoS        = 8906.26;          % [Kg/m^3]
    Tinf        = FarTemp(k);       % [K]
    Vinf        = 50;               % [m/s]
    Pinf        = 9e4;              % [Pa]
    Tfreezing   = 273.15;           % [K]
    Lambda_e    = 2.26e6;           % [J/Kg]
    Lambda_f    = 3.344e5;          % [J/Kg]
    Lambda_s    = 2.83e3;           % [J/Kg]
    cpAir       = 1014.5;           % [J/(Kg*K)]
    cpWater     = 4220;             % [J/(Kg*K)]
    cpS         = 385.19;           % [J/(Kg*K)]
    kWater      = 0.571;            % [W/(m*K)]
    kIce        = 2.18;             % [W/(m*K)]
    if (IPS && t_IPS == 0)
        q_IPS   = IPS_Heat_Flux;    % [W/m^2]
    else
        q_IPS   = 0;
    end


    t0          = 0;                % [s]
    tEnd        = 54;              % [s]
    dt          = 0.1;              % [s]

    % I dati del substrato sono quelli dell'heating element della tesi
    %% Initialize Variables

    q_aero      = 0;
    q_kin       = 0;
    q_conv      = 0;
    q_drop      = 0;
    q_evap      = 0;
    q_sub       = 0;
    h0          = 0;
    B0          = 0;
    c0          = 0;
    h           = h0;
    B           = B0;
    c           = c0;
    Twater      = Twall;
    TwaterMelt  = Tfreezing;
    Tdrop       = Tinf;
    Tice        = Tfreezing;
    IceType     = 0; % 0-AS, 1-AWS, 2-AIS, 3-AWIS, 4-AIWS, 5-AWIWS
    IceRate     = 0;

    
    % Compute average quantities on the substratum, will be usefull also
    % when including conduction in the substratum 
    
    H_H     = 0.03*10^-3;
    H_Er    = 0.20*10^-3;
    H_El    = 0.26*10^-3;
    H_Fg    = 0.89*10^-3;
    H_Sil   = 3.43*10^-3;

    rCp_H   = 8906.26 * 385.19;
    rCp_Er  = 8025.25 * 502.42;
    rCp_El  = 1384 * 1256.04;
    rCp_Fg  = 1794.07 * 1570.05;
    rCp_Sil = 648.25 * 125.6;
    
    k_H     = 41.02; %[W/(m*K)] 
    k_Er    = 16.27;
    k_El    = 0.256;
    k_Fg    = 0.294;
    k_Sil   = 0.121;
    
    ks       = (k_H * H_H + k_Er * H_Er + 2 * k_El * H_El + k_Fg * H_Fg + k_Sil *  H_Sil) / (H_H + H_Er + 2*H_El + H_Fg + H_Sil); 
    rhoCpS  = (rCp_H * H_H + rCp_Er * H_Er + 2 * rCp_El * H_El + rCp_Fg * H_Fg + rCp_Sil *  H_Sil) / (H_H + H_Er + 2*H_El + H_Fg + H_Sil); 
    H       = H_H + H_Er + 2*H_El + H_Fg + H_Sil;



    %% Start simulation
    maxVal = 0; % Variable used in plot height to define the upper limit of plot
    j   = 1; % Variable used in plot height to define the current number of the plot 
    ind = 0; % Variable used as index to store the height of the ice, liquid film and static liquid film
    tol = 1e-8; % Defined tolerance for convergence 
    itmax = 12; % Defined max iterations
    
    % Initialize variables to compute the error
    Bold = B0;
    hOld = h0;
    cold = c0;
    err = 10;
    it = 0;
    
    
    Temperatura_parete = [];
    Temperatura_acqua0 = [];
    Heat_Flux = [];
%% Initial condition, water impinge on the clean surface: AS
    while (err > tol && it < itmax)
        IceType = AS();

        if Bold == 0
            denB = 1;
        end
        if hOld == 0
            denh = 1;
        end
        if cold == 0
            denc = 1;
        end

        err = max([(B-Bold)/(denB),(h-hOld)/(denh), (c-cold)/(denc)]);
        
        Bold = B;
        hOld = h;
        cold = c;

        it = it + 1;
        
        if (PrintErrorIt)
            fprintf('\n Err \t = %1.10f \t It \t = %d \n ', err, it);
        end
        
        if (PrintValuesIt)
            fprintf('\n Twall \t = %3.2f \n B \t = %1.10f \n h \t = %1.10f \n c \t = %1.10f \n dh \t = %1.10f \n dB \t = %1.10f \n dc \t = %1.10f \n \n ', Twall, B,h,c, h-h0, B - B0, c-c0);
        end
    end

        
        if (PrintErrorEnd)
            fprintf('\n Err \t = %1.10f \t It \t = %d \n ', err, it);
        end
        
        if (PrintValuesEnd)
            fprintf('\n Twall \t = %3.2f \n B \t = %1.10f \n h \t = %1.10f \n c \t = %1.10f \n dh \t = %1.10f \n dB \t = %1.10f \n dc \t = %1.10f \n \n ', Twall, B,h,c, h-h0, B - B0, c-c0);
        end

    
    
    % Save solution first time step
    ind = ind + 1;
    B_tot(ind) = B;
    h_tot(ind) = h;
    c_tot(ind) = c;
    M_SUB(ind) = m_sub;
    M_IMP(ind) = m_imp;
    M_EVAP(ind) = m_evap;
    M_ICE(ind) = m_ice + m_melt; %rhoIceGlaze*(B - B0)/dt; 
    M_WATER(ind) = m_water; %rhoWater*(h-h0)/dt
    M_MELT(ind) = m_waterMelt;
%     if c-c0 >= 0
%         M_MELT(ind) = m_waterMelt %rhoWater*(c-c0)/dt;
%     else
%         M_MELT(ind) = 0;
%         M_ICE(ind) = 0;
%     end
    
    % Update initial solution for iterations and time step
    B0 = B;
    h0 = h;
    c0 = c;
    
    t = dt + t0;

    if (PlotTempHeights)
        PlotNum(j) = plotHeights(t, IceType,j);
        j = j + 1;
    end


    %Update wall temperature increased caused by IPS
     %Twall = Twall + q_IPS * dt / (rhoCpS * H);
     
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       N = 50;  %50 se usiamo la Twall_cond, 70 se usiamo Twall_cond_strat
       Temp_wall = Twall*ones(1,N+2);
      [Temp_wall, Q] = Twall_cond_test( H, rhoCpS, ks, q_IPS, Temp_wall, IceType,Twater);
      %[Temp_wall] = Twall_cond_strat(q_IPS, Temp_wall, IceType, Twall);
      Twall = Temp_wall(end); 
     
      Temperatura_parete = [Temperatura_parete Twall];
      Temperatura_acqua0 = [Temperatura_acqua0 Twall];
      Heat_Flux = [Heat_Flux Q];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    if B + h + c > maxVal
            maxVal = B + h + c;
    end
    
    
%% Second time step, start ice type can be: AS, AWS, AIS, AWIS 
    err = 10;
    it  = 0;

    fprintf('\n t = %f Tipo = %d \n \n ', t, IceType);


    while(err > tol && it < itmax)
    if IceType == 0
        hAir =  1000;
        IceType = AS();
    elseif IceType == 1
        hAir =  500;
        IceType = AWS();
    elseif IceType == 2
        hAir =  1000;
        IceType = AIS();
    elseif IceType == 3
        hAir =  500;
        IceType = AWIS();
    end

        if Bold == 0
            denB = 1;
        end
        if hOld == 0
            denh = 1;
        end
        if cold == 0
            denc = 1;
        end

        err = max([(B-Bold)/(denB),(h-hOld)/(denh), (c-cold)/(denc)]);

        Bold = B;
        hOld = h;
        cold = c;
        it = it + 1;
        
        if (PrintErrorIt)
            fprintf('\n Err \t = %1.10f \t It \t = %d \n ', err, it);
        end
        
        if (PrintValuesIt)
            fprintf('\n Twall \t = %3.2f \n B \t = %1.10f \n h \t = %1.10f \n c \t = %1.10f \n dh \t = %1.10f \n dB \t = %1.10f \n dc \t = %1.10f \n \n ', Twall, B,h,c, h-h0, B - B0, c-c0);
        end
    end

        
        if (PrintErrorEnd)
            fprintf('\n Err \t = %1.10f \t It \t = %d \n ', err, it);
        end
        
        if (PrintValuesEnd)
            fprintf('\n Twall \t = %3.2f \n B \t = %1.10f \n h \t = %1.10f \n c \t = %1.10f \n dh \t = %1.10f \n dB \t = %1.10f \n dc \t = %1.10f \n \n ', Twall, B,h,c, h-h0, B - B0, c-c0);
        end

    % Save solution second time step
    ind = ind + 1;
    B_tot(ind) = B;
    h_tot(ind) = h;
    c_tot(ind) = c;
    
    M_SUB(ind) = m_sub;
    M_IMP(ind) = m_imp;
    M_EVAP(ind) = m_evap;
    M_ICE(ind) = m_ice + m_melt; %rhoIceGlaze*(B - B0)/dt; 
    M_WATER(ind) = m_water; %rhoWater*(h-h0)/dt
    M_MELT(ind) = m_waterMelt;
%     if c-c0 >= 0
%         M_MELT(ind) = rhoWater*(c-c0)/dt;
%     else
%         M_MELT(ind) = 0;
%         M_ICE(ind) = 0;
%     end
    
    
    % Update initial solution for iterations and time step
    B0 = B;
    h0 = h;
    c0 = c;
    
    t = dt + t;

    if (PlotTempHeights)
        PlotNum(j) = plotHeights(t, IceType,j);
        j = j + 1;
    end


    % Update wall temperature increased caused by IPS
     %Twall = Twall + q_IPS * dt / (rhoCpS * H);
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [Temp_wall, Q] = Twall_cond_test( H, rhoCpS, ks, q_IPS, Temp_wall, IceType);
      %[Temp_wall] = Twall_cond_strat(q_IPS, Temp_wall, IceType, Twall);
      Twall = Temp_wall(end); 
      
     Temperatura_parete = [Temperatura_parete Twall];
     Heat_Flux = [Heat_Flux Q];
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if B + h + c > maxVal
        maxVal = B + h + c;
    end

%% Following time steps, start ice type can be: AS, AWS, AIS, AWIS, AIWS, AWIWS
    while t < tEnd
        
     err = 10;
     it  = 0;
     

     while (err > tol && it < itmax)
         
        if IceType == 0
            hAir =  1000;
            IceType = AS_mod();
        elseif IceType == 1 
            hAir =  500;
            IceType = AWS();
        elseif IceType == 2
            hAir =  1000;
            IceType = AIS();
        elseif IceType == 3
            hAir =  500;
            IceType = AWIS();
        elseif IceType == 4
            hAir =  1000;
            IceType = AIWS();
        elseif IceType == 5
            hAir =  500;
            IceType = AWIWS();
        end

        if Bold == 0
            denB = 1;
        end
        if hOld == 0
            denh = 1;
        end
        if cold == 0
            denc = 1;
        end

        err = max([(B-Bold)/(denB),(h-hOld)/(denh), (c-cold)/(denc)]);
        
        it = it + 1;

        Bold = B;
        hOld = h;
        cold = c;
        
        if (PrintErrorIt)
            fprintf('\n Err \t = %1.10f \t It \t = %d \n ', err, it);
        end
        
        if (PrintValuesIt)
            fprintf('\n Twall \t = %3.2f \n B \t = %1.10f \n h \t = %1.10f \n c \t = %1.10f \n dh \t = %1.10f \n dB \t = %1.10f \n dc \t = %1.10f \n \n ', Twall, B,h,c, h-h0, B - B0, c-c0);
        end
    end

        
        if (PrintErrorEnd)
            fprintf('\n Err \t = %1.10f \t It \t = %d \n ', err, it);
        end
        
        if (PrintValuesEnd)
            fprintf('\n Twall \t = %3.2f \n B \t = %1.10f \n h \t = %1.10f \n c \t = %1.10f \n dh \t = %1.10f \n dB \t = %1.10f \n dc \t = %1.10f \n \n ', Twall, B,h,c, h-h0, B - B0, c-c0);
        end 
        
    ind = ind + 1;
    B_tot(ind) = B;
    h_tot(ind) = h;
    c_tot(ind) = c;
    
    M_SUB(ind) = m_sub;
    M_IMP(ind) = m_imp;
    M_EVAP(ind) = m_evap;
        M_ICE(ind) = m_ice + m_melt; %rhoIceGlaze*(B - B0)/dt; 
    M_WATER(ind) = m_water; %rhoWater*(h-h0)/dt;
    M_MELT(ind) = m_waterMelt;
%     if c-c0 >= 0
%         M_MELT(ind) = rhoWater*(c-c0)/dt;
%     else
%         M_MELT(ind) = 0;
%         M_ICE(ind) = 0;
%     end
    
    % Update initial solution for iterations and time step
    B0 = B;
    h0 = h;
    c0 = c;
    
    t = dt + t;

    if (PlotTempHeights)
        if (mod(t,1) < 1e-12)   || mod(t,1) > 0.999
            PlotNum(j) = plotHeights(t, IceType,j);
            j = j + 1;
        end
    end
    if (IPS)
        if t > t_IPS
            q_IPS = IPS_Heat_Flux;
        end
    end

    %Update wall temperature increased caused by IPS
    % Twall = Twall + q_IPS * dt / (rhoCpS * H);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [Temp_wall, Q] = Twall_cond_test( H, rhoCpS, ks, q_IPS, Temp_wall, IceType);
      %[Temp_wall] = Twall_cond_strat(q_IPS, Temp_wall, IceType,Twall);
      Twall = Temp_wall(end); 
      
      Temperatura_parete = [Temperatura_parete Twall];
      Heat_Flux = [Heat_Flux Q];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
        if B + h + c > maxVal
            maxVal = B + h + c;
        end

    end
    
    if (PlotTempHeights)
        for i = 1:j-1
            PlotNum(i).YLim = [-3e-5 maxVal*1.2]; 
        end
    end
    
    if (PlotTotal) 
        fS = figure('Position', get(0, 'Screensize'));
        subplot(2,1,1)
        hold on
        if length([t0:dt:tEnd]) == length(B_tot)
            plot ([t0:dt:tEnd], [B_tot], 'LineWidth',l, 'DisplayName',['B'])
            plot ([t0:dt:tEnd],[h_tot], 'LineWidth',l, 'DisplayName',[' h'])
            plot ([t0:dt:tEnd],[c_tot ], 'LineWidth',l, 'DisplayName',['c '])
            plot ([t0:dt:tEnd],[B_tot+c_tot+h_tot ], '--','LineWidth',l, 'DisplayName',['B + h + c'])
        else
            plot ([t0:dt:tEnd], [0, B_tot], 'LineWidth',l, 'DisplayName',['B'])
            plot ([t0:dt:tEnd],[0, h_tot], 'LineWidth',l, 'DisplayName',[' h'])
            plot ([t0:dt:tEnd],[0, c_tot ], 'LineWidth',l, 'DisplayName',['c '])
            plot ([t0:dt:tEnd],[0, B_tot+c_tot+h_tot ], '--','LineWidth',l, 'DisplayName',['B + h + c'])
        end

%         xlim([0 10])
%         ylim([-0.035 0.025])
        set(gca,'FontSize',fM);
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$t [s]$', 'Interpreter','latex')
%         ylabel('$ [m]$', 'Interpreter','latex')
        hl = legend('show');
        legend boxoff
        set(hl,'Location','westoutside', 'Interpreter','latex','FontSize',fM)

        subplot(2,1,2)
        hold on
        if length([t0:dt:tEnd]) == length(B_tot)
            plot ([t0:dt:tEnd], [M_IMP], 'LineWidth',l, 'DisplayName',['$m_{imp}$'])
            plot ([t0:dt:tEnd],[M_EVAP], 'LineWidth',l, 'DisplayName',['$m_{evap}$'])
            plot ([t0:dt:tEnd],[M_SUB], 'LineWidth',l, 'DisplayName',['$m_{sub}$'])
            plot ([t0:dt:tEnd],[M_ICE], 'LineWidth',l, 'DisplayName',['$m_{ice}$'])
            plot ([t0:dt:tEnd],[M_WATER], 'LineWidth',l, 'DisplayName',['$m_{water}$'])
            plot ([t0:dt:tEnd],[M_MELT], 'LineWidth',l, 'DisplayName',['$m_{melt}$'])
            plot ([t0:dt:tEnd],[M_MELT + M_EVAP + M_ICE + M_WATER + M_SUB], '--', 'LineWidth',l, 'DisplayName',['$\Sigma m$'])
        else
           plot ([t0:dt:tEnd], [0, M_IMP], 'LineWidth',l, 'DisplayName',['$m_{imp}$'])
            plot ([t0:dt:tEnd],[0, M_EVAP], 'LineWidth',l, 'DisplayName',['$m_{evap}$'])
            plot ([t0:dt:tEnd],[0, M_SUB], 'LineWidth',l, 'DisplayName',['$m_{sub}$'])
            plot ([t0:dt:tEnd],[0, M_ICE], 'LineWidth',l, 'DisplayName',['$m_{ice}$'])
            plot ([t0:dt:tEnd],[0, M_WATER], 'LineWidth',l, 'DisplayName',['$m_{water}$'])
            plot ([t0:dt:tEnd],[0, M_MELT], 'LineWidth',l, 'DisplayName',['$m_{melt}$'])
            plot ([t0:dt:tEnd],[0, M_MELT + M_EVAP + M_ICE + M_WATER + M_SUB],'--', 'LineWidth',l, 'DisplayName',['$\Sigma m$'])
        end

        set(gca,'FontSize',fM);
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$t [s]$', 'Interpreter','latex')
        hl = legend('show');
        legend boxoff
        set(hl,'Location','westoutside', 'Interpreter','latex','FontSize',fM)

        if save
            saveas(fS, [pathSave 'LessWater_IPS_10_' num2str(Twall) '_' num2str(Tinf) '.eps'], 'epsc');
        end
    end



    end
end

figure
plot(t0:dt:tEnd,Temperatura_parete)
hold on
plot(t0:dt:tEnd,Temperatura_acqua0)
xlabel('Time [s]'); ylabel('Temperature [K]');
title('Twall')
figure 
plot(t0:dt:tEnd,Heat_Flux)