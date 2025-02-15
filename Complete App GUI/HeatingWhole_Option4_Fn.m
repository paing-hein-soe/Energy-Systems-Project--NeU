function [T2A, T2B, T4A, P2A, P2B, P4A, ...
    ...
    Energy_w_rate_in_hp, Energy_q_rate_in_hp_free, Energy_q_rate_out_hp, ...
    Energy_w_rate_out_he, Energy_q_rate_in_he, Energy_q_rate_out_he, ...
    Energy_q_rate_loss_furnace, Energy_q_rate_in_furnace, Energy_q_rate_out_furnace,...
    Energy_q_rate_out_regHE, Exg_f_rate_out_regHE,... %% Additional for Option 4
    ...
    Energy_q_rate_in_heatingwhole, Energy_q_rate_out_heatingwhole,...
    Energy_w_rate_in_heatingwhole,Energy_w_rate_out_heatingwhole,...
    Energy_q_rate_in_heatingwhole_free, Energy_q_rate_loss_heatingwhole,...
    Energy_rate_in_heatingwhole, Energy_utilized_rate_heatingwhole,...
    ...
    Exg_f_rate_out_3, Exg_f_rate_in_4,...
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
    eff_energy_heatingwhole, eff_exergy_heatingwhole] = HeatingWhole_Option4_Fn(T0,T2,T3,T4,T_flame,...
    P0, P2, P3, P4, m_dot, cp, R ,cop_hp_factor, eta_he_factor, eta_furnace)
%% Whole System of Furnace Heating - Heating Air + Furnace

%% Energy Analysis
%% Finding T2 A (Before first heating from HE) and T2B (Before heating from HP)
% T2A using regenerative heat exchangeer
eta_regHE = 0.8;
T2A =T2 +  eta_regHE * (T4- T2); % Heat Exchanger before Heat Engine
dT = T2A-T2; T4A = T4 - dT;
%
cop_hp_ideal = T3/(T3-T0);
cop_hp_actual = cop_hp_ideal * cop_hp_factor; % cop_hp_factor of ideal
                                                                            
eta_he_ideal = (T_flame-T0)/T_flame;
eta_he_actual= eta_he_ideal * eta_he_factor; % cop_he_factor of ideal

F=cop_hp_actual * 1/(1-eta_he_actual)* eta_he_actual;

T2B = (F*T2A +T3)/(1+F);
%
P2A = P2; P2B=P2; P4A = P4; % No Pressure change

if T4<= T2
    error('Discharge Temperature T4 (at the exit of Hospital) is less than or equal to T2 (at Fan Exit)')
end
%% Calculations for Regenerative Heat Exchanger
Energy_q_rate_in_heatingair_regHE = m_dot*cp*(T2A-T2);
Energy_q_rate_out_regHE = m_dot*cp*(T4-T4A); % The Same with Energy_q_rate_in_heatingair_regHE

%% Calculations for Output Heat Rate Required for HP and HE 
Energy_q_rate_in_heatingair_hp  = m_dot*cp*(T3-T2B);
Energy_q_rate_in_heatingair_he = m_dot*cp*(T2B-T2A);
%% 
Energy_q_rate_out_hp = Energy_q_rate_in_heatingair_hp; %Linked to Heating Air Sytem
Energy_w_rate_in_hp = Energy_q_rate_out_hp/cop_hp_actual;
Energy_q_rate_in_hp_free = Energy_q_rate_out_hp - Energy_w_rate_in_hp; % FREE source

Energy_q_rate_out_he = Energy_q_rate_in_heatingair_he; %Linked to Heating Air Sytem
Energy_w_rate_out_he = Energy_w_rate_in_hp; %Linked to HP
Energy_q_rate_in_he = Energy_w_rate_out_he/eta_he_actual;
Exg_q_rate_loss_he = 0;
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
Energy_q_rate_loss_heatingwhole = Energy_q_rate_loss_furnace;

Energy_rate_in_heatingwhole = Energy_q_rate_in_heatingwhole + Energy_w_rate_in_heatingwhole;  %Does not include free heat from HP
Energy_utilized_rate_heatingwhole = m_dot * cp* (T3-T2); % For example - Take T3 -T2

%% Exergy Analysis
% Heating Air

s3_s2 = cp * log(T3 / T2) - R * log(P3 / P2); % Entropy change (J/kg·K)
s4A_s4 = cp * log(T4A / T4) - R * log(P4A / P4); % Entropy change (J/kg·K)

s2_s0 = cp * log(T2 / T0) - R * log(P2/ P0); % Entropy change (J/kg·K)
s3_s0 = cp * log(T3 / T0) - R * log(P3 / P0); % Entropy change (J/kg·K)
s4_s0 = cp * log(T4 / T0) - R * log(P4/ P0); % Entropy change (J/kg·K)
s4A_s0 = cp * log(T4A / T0) - R * log(P4A/ P0); % Entropy change (J/kg·K)

exg_f_in_2 = cp * (T2 - T0) - T0 * s2_s0; %exg_f_out_fan; 
exg_f_out_3 = cp * (T3 - T0) - T0 * s3_s0; %exg_f_out_fan; 
exg_f_in_4 = cp * (T4 - T0) - T0 * s4_s0; %exg_f_out_fan; 
exg_f_out_4A = cp * (T4A - T0) - T0 * s4A_s0; %exg_f_out_fan; 

Exg_f_rate_in_2A = m_dot * exg_f_in_2;
Exg_f_rate_in_4 = m_dot * exg_f_in_4;
Exg_f_rate_out_3 = m_dot * exg_f_out_3;
Exg_f_rate_out_4A = m_dot * exg_f_out_4A;

Exg_f_rate_out_regHE = Exg_f_rate_out_4A;
% HP
Exg_q_rate_in_hp_free = (1-T0/T0)*Energy_q_rate_in_hp_free; % T0 - Assumptions 
% HE %Tb_HE = T0; % Assumptions T0 %Exg_q_rate_loss_he = (1-T0/Tb_HE)*Energy_q_rate_loss_he; 
% Furnace
Exg_q_rate_in_furnace = (1-T0/T_flame)*Energy_q_rate_in_furnace;
Exg_q_rate_loss_furnace = (1-T0/T_flame)*Energy_q_rate_loss_furnace ; % T_flame Assumptions

%% Exergy Destroyed from Entropy Balance 
sigma_dot_p_heatingwhole =-(Energy_q_rate_in_hp_free/T0...
    -Energy_q_rate_loss_furnace/T_flame +Energy_q_rate_in_furnace/T_flame)+ m_dot * s3_s2 + m_dot *s4A_s4;

Exg_d_rate_heatingwhole = T0*sigma_dot_p_heatingwhole;

%% Exergy Destroyed from Exergy Balance
exg_f_in_heatingwhole = exg_f_in_2 + exg_f_in_4; %Summation
exg_f_out_heatingwhole = exg_f_out_3 + exg_f_out_4A; % Summation

Exg_f_rate_in_heatingwhole = m_dot * exg_f_in_heatingwhole;
Exg_f_rate_out_heatingwhole = m_dot * exg_f_out_heatingwhole;

Exg_q_rate_in_heatingwhole = Exg_q_rate_in_furnace ;
Exg_q_rate_out_heatingwhole = 0; 

Exg_w_rate_in_heatingwhole = 0;
Exg_w_rate_out_heatingwhole = 0;

Exg_q_rate_in_heatingwhole_free = Exg_q_rate_in_hp_free;
Exg_q_rate_loss_heatingwhole = Exg_q_rate_loss_furnace ;

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
fprintf('\nThe Whole Heating System - Furnace + (RegHeatExchanger+HE+HP) + Heating Air- Results:\n');

fprintf('Exergy Destroyed Rate in the whole Heating System (W): %.2f\n', Exg_d_rate_heatingwhole);
fprintf('Exergy Destroyed Rate (Check- Exergy Balance) in the whole Heating System (W): %.2f\n', Exg_d_rate_heatingwhole_exgbalance);

end