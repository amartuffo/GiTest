function IceType = AS_mod()

global Beta LWC r hAir hWater Twall rhoIceRime rhoIceGlaze rhoWater Tinf Vinf Tfreezing Lambda_e Lambda_f Lambda_s cpAir cpWater kWater kIce dt;
global h0 B0 h B c Twater Tdrop Tice
global m_imp m_evap m_sub m_water m_ice m_melt m_waterMelt BC

% Evaluates heat fluxes 
q_evap  = EvapHeatFlux(Twall, Tinf);
q_aero  = 0.5 * r * hAir * Vinf^2 / cpAir;
q_kin   = 0.5 * LWC * Beta * Vinf^3;
q_conv  = hAir * (Twall - Tinf);
q_drop  = LWC * Beta * Vinf * cpWater * (Twall - Tdrop);
q_sub   = 0;

m_imp   = Beta * LWC * Vinf;
m_evap  = q_evap / Lambda_e;
m_sub   = 0;
m_ice   = 0;
m_water = 0;
m_melt  = 0;
m_waterMelt = 0;

% Compute water film height as in the AWS case, if h<0 we are still in AS
% otherwise either water or ice is forming.
WaterRate = AWSWaterRate(m_imp, m_evap, rhoWater);
h         = WaterRate * dt + h0;

if (h < 0 )
    % All the impinging water is evaporating, TWater will be = Twall and
    % the if statements will end on the fully evaporative
    h = 0;
    m_evap = m_imp;
    q_evap  = m_evap * Lambda_e;
end

Twater  = TWaterAWS(q_aero, q_kin, h/2 ,BC);

if Twater < Tfreezing 
    fprintf('\n Ice: ');
    % Ice is forming
    IceRate   = GlazeIceRate (q_aero, q_kin, q_conv, q_drop, q_evap, rhoIceGlaze, Lambda_f, B, kIce);
    m_ice     = rhoIceGlaze * IceRate;
    B         = IceRate * dt + B0;
    WaterRate = GlazeWaterRate(m_imp, m_evap,m_ice, rhoWater);
    h         = WaterRate * dt + h0;
    if h < 0 || IceRate < 0
        % Rime ice is forming
        fprintf(' --> Rime Ice \n');
        q_sub   = SubHeatFlux(Twall, Tinf);
        m_sub   = q_sub / Lambda_s;
        IceRate = RimeIceRate(m_imp, m_sub, rhoIceRime);
        m_ice   = rhoIceRime * IceRate;
        q_lat   = m_ice * Lambda_f;
        B       = IceRate * dt + B0;
        h       = 0;
        c       = 0;
        Tice    = TIceRime (q_aero, q_kin, q_lat, B ,BC);
        Twater  = Tfreezing;
        m_evap  = 0;
        m_water = 0;
        m_melt  = 0;
        IceType = 2;
    else
        fprintf(' --> Glaze Ice \n');
        Twater   = TWaterGlaze (q_aero, q_kin, h);
        Tice     = TIceGlaze(Twall, Tfreezing, B, B/2 ,BC);
        m_sub    = 0;
        m_melt   = 0;
        m_water  = rhoWater * WaterRate;
        IceType  = 3;
    end
else
    fprintf('\n No Ice: ');
    if h > 1e-5
        fprintf(' --> Running Wet \n');
        WaterRate = AWSWaterRate(m_imp, m_evap, rhoWater);
        h         = WaterRate * dt + h0;
        B         = 0;
        c         = 0;
        m_sub   = 0;
        m_ice   = 0;
        m_melt  = 0;
        m_water = rhoWater * WaterRate;
        Twater  = TWaterAWS (q_aero, q_kin, h ,BC);
        Tice    = Tfreezing;
        IceType   = 1;
    else
        fprintf(' --> Fully Evaporative \n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WaterRate = AWSWaterRate(m_imp, m_evap, rhoWater);
        m_evap  = m_evap + rhoWater*WaterRate;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        B       = 0;
        c       = 0;
        h       = 0;
        m_sub   = 0;
        m_ice   = 0;
        m_water = 0;
        m_melt  = 0;
        Twater  = Twall;
        Tice    = Tfreezing;
        IceType = 0;
    end
end

end