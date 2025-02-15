function [...
    Energy_w_rate_in_hp, Energy_q_rate_in_hp_free, Energy_q_rate_out_hp, ...
    Energy_w_rate_out_he, Energy_q_rate_in_he, Energy_q_rate_loss_he, ...
    Energy_q_rate_loss_furnace, Energy_q_rate_in_furnace, Energy_q_rate_out_furnace, ...
    ...
    Energy_q_rate_in_heatingwhole, Energy_q_rate_out_heatingwhole,...
    Energy_w_rate_in_heatingwhole,Energy_w_rate_out_heatingwhole,...
    Energy_q_rate_in_heatingwhole_free, Energy_q_rate_loss_heatingwhole,...
    Energy_rate_in_heatingwhole, Energy_utilized_rate_heatingwhole,...
    ...
    exg_f_in_heatingwhole, exg_f_out_heatingwhole,...
    Exg_f_rate_in_heatingwhole, Exg_f_rate_out_heatingwhole,...
    Exg_q_rate_in_heatingwhole, Exg_q_rate_out_heatingwhole,...
    Exg_w_rate_in_heatingwhole, Exg_w_rate_out_heatingwhole,...
    Exg_q_rate_in_heatingwhole_free, Exg_q_rate_loss_heatingwhole,...
    Exg_rate_in_heatingwhole, Exg_utilized_rate_heatingwhole,...
    Exg_q_rate_loss_furnace,Exg_q_rate_loss_he,Exg_q_rate_in_hp_free,... % Additional for Option 2
    ...
    sigma_dot_p_heatingwhole, Exg_d_rate_heatingwhole,Exg_d_rate_heatingwhole_exgbalance,...
    eff_energy_heatingwhole, eff_exergy_heatingwhole] = HeatingWhole_Option2_Fn(T0,T2,T3,T_flame,...
    P0, P2, P3, m_dot, cp, R ,cop_hp_factor, eta_he_factor, eta_furnace, Q_dot_heatingair, exg_f_out_fan )
%% Whole System of Furnace Heating - Heating Air + Furnace

%% Energy Analysis
% Heating Air 
Energy_q_rate_in_heatingair = Q_dot_heatingair;

% Heat Pump
cop_hp_ideal = T3/(T3-T0);
cop_hp_actual = cop_hp_ideal * cop_hp_factor; % cop_hp_factor of ideal

Energy_q_rate_out_hp = Energy_q_rate_in_heatingair; %Linked to Heating Air Sytem
Energy_w_rate_in_hp = Energy_q_rate_out_hp/cop_hp_actual;
Energy_q_rate_in_hp_free = Energy_q_rate_out_hp - Energy_w_rate_in_hp; % FREE source

% Heat Engine
eta_he_ideal = (T_flame-T0)/T_flame;
eta_he_actual= eta_he_ideal * eta_he_factor; % cop_he_factor of ideal

Energy_w_rate_out_he = Energy_w_rate_in_hp;
Energy_q_rate_in_he = Energy_w_rate_out_he/eta_he_actual;
Energy_q_rate_loss_he =Energy_q_rate_in_he-Energy_w_rate_out_he;

% Furnace
Energy_q_rate_out_furnace = Energy_q_rate_in_he; %Linked to he
Energy_q_rate_in_furnace = Energy_q_rate_out_furnace/ eta_furnace; % Heat input from fuel (W)
Energy_q_rate_loss_furnace = Energy_q_rate_in_furnace-Energy_q_rate_out_furnace;

% Heating whole << ALL Properties >>
Energy_q_rate_in_heatingwhole = Energy_q_rate_in_furnace; %Does not include free heat input in HP
Energy_q_rate_out_heatingwhole = 0;
Energy_w_rate_in_heatingwhole = 0;
Energy_w_rate_out_heatingwhole = 0;

Energy_q_rate_in_heatingwhole_free = Energy_q_rate_in_hp_free; 
Energy_q_rate_loss_heatingwhole = Energy_q_rate_loss_he+Energy_q_rate_loss_furnace;

Energy_rate_in_heatingwhole = Energy_q_rate_in_heatingwhole + Energy_w_rate_in_heatingwhole;  %Does not include free heat from HP
Energy_utilized_rate_heatingwhole = m_dot * cp* (T3-T2);

%% Exergy Analysis
% Heating Air
s3_s2 = cp * log(T3 / T2) - R * log(P3 / P2); % Entropy change (J/kg·K)
s3_s0 = cp * log(T3 / T0) - R * log(P3 / P0); % Entropy change (J/kg·K)

exg_f_in_heatingair = exg_f_out_fan; % Specific exergy at Heating by Furnace heating inlet (J/kg)
exg_f_out_heatingair= cp * (T3 - T0) - T0 * s3_s0; % Specific exergy at outlet (J/kg)
Exg_f_rate_in_heatingair = m_dot * exg_f_in_heatingair;
Exg_f_rate_out_heatingair = m_dot * exg_f_out_heatingair;
% HP
Exg_q_rate_in_hp_free = (1-T0/T0)*Energy_q_rate_in_hp_free; % T0 - Assumptions 
% HE
Tb_HE = T0; % Assumptions T0
Exg_q_rate_loss_he = (1-T0/Tb_HE)*Energy_q_rate_loss_he; 
% Furnace
Exg_q_rate_in_furnace = (1-T0/T_flame)*Energy_q_rate_in_furnace;
Exg_q_rate_loss_furnace = (1-T0/T_flame)*Energy_q_rate_loss_furnace ; % T_flame Assumptions

%% Exergy Destroyed from Entropy Balance 
sigma_dot_p_heatingwhole =-(Energy_q_rate_in_hp_free/T0 - Energy_q_rate_loss_he/Tb_HE ...
    -Energy_q_rate_loss_furnace/T_flame +Energy_q_rate_in_furnace/T_flame)+ m_dot * s3_s2;

Exg_d_rate_heatingwhole = T0*sigma_dot_p_heatingwhole;

%% Exergy Destroyed from Exergy Balance
exg_f_in_heatingwhole = exg_f_in_heatingair;
exg_f_out_heatingwhole = exg_f_out_heatingair;

Exg_f_rate_in_heatingwhole = Exg_f_rate_in_heatingair;
Exg_f_rate_out_heatingwhole = Exg_f_rate_out_heatingair;

Exg_q_rate_in_heatingwhole = Exg_q_rate_in_furnace ;
Exg_q_rate_out_heatingwhole = 0; 

Exg_w_rate_in_heatingwhole = 0;
Exg_w_rate_out_heatingwhole = 0;

Exg_q_rate_in_heatingwhole_free = Exg_q_rate_in_hp_free;
Exg_q_rate_loss_heatingwhole = Exg_q_rate_loss_furnace + Exg_q_rate_loss_he;

Exg_rate_in_heatingwhole = Exg_q_rate_in_heatingwhole+Exg_w_rate_in_heatingwhole;
Exg_utilized_rate_heatingwhole=Exg_f_rate_out_heatingwhole-Exg_f_rate_in_heatingwhole;

Exg_d_rate_heatingwhole_exgbalance= (Exg_q_rate_in_heatingwhole-Exg_q_rate_out_heatingwhole)...
   + (Exg_q_rate_in_heatingwhole_free -Exg_q_rate_loss_heatingwhole)+...
   (Exg_w_rate_in_heatingwhole-Exg_w_rate_out_heatingwhole)+...
   (Exg_f_rate_in_heatingwhole -Exg_f_rate_out_heatingwhole);


%% Efficiency
eff_energy_heatingwhole= Energy_utilized_rate_heatingwhole/Energy_rate_in_heatingwhole;
eff_exergy_heatingwhole= Exg_utilized_rate_heatingwhole/Exg_rate_in_heatingwhole;

%% Energy/Exergy Results -Heating
fprintf('\nThe Whole Heating System - Furnace + (HE+HP) + Heating Air- Results:\n');

fprintf('Exergy Destroyed Rate in the whole Heating System (W): %.2f\n', Exg_d_rate_heatingwhole);
fprintf('Exergy Destroyed Rate (Check- Exergy Balance) in the whole Heating System (W): %.2f\n', Exg_d_rate_heatingwhole_exgbalance);

end