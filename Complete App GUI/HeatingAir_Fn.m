function [Energy_w_rate_in_heatingair, Energy_q_rate_in_heatingair, Energy_rate_in_heatingair, ...
    Energy_w_rate_out_heatingair, Energy_q_rate_out_heatingair, Energy_rate_out_heatingair, Energy_q_rate_loss_heatingair,...
    ...
    Exg_w_rate_in_heatingair, Exg_q_rate_in_heatingair, Exg_rate_in_heatingair, ...
    Exg_w_rate_out_heatingair, Exg_q_rate_out_heatingair, Exg_rate_out_heatingair, Exg_q_rate_loss_heatingair, ...
    ...
    exg_f_in_heatingair, exg_f_out_heatingair, Exg_f_rate_in_heatingair, Exg_f_rate_out_heatingair, ...
    sigma_dot_p_heatingair, Exg_d_rate_heatingair, Exg_d_rate_heatingair_exgbalance, ...
    ...
    Energy_utilized_rate_heatingair, Exg_utilized_rate_heatingair, eff_exergy_heatingair, eff_energy_heatingair]=HeatingAir_Fn...
    (m_dot, cp, R, T0, T2, T3, P0, P2, P3, Q_dot_heatingair, T_flame, exg_f_out_fan )
% Heating Air System (Air Flow from state 2 to 3) -- Primary - 

%% Energy Analysis for Heating Air SYstem

Energy_q_rate_in_heatingair = Q_dot_heatingair;
Energy_w_rate_in_heatingair = 0;
Energy_rate_in_heatingair= Energy_q_rate_in_heatingair+Energy_w_rate_in_heatingair;

Energy_q_rate_out_heatingair= 0;
Energy_w_rate_out_heatingair = 0;
Energy_rate_out_heatingair= Energy_q_rate_out_heatingair + Energy_w_rate_out_heatingair;

Energy_utilized_rate_heatingair = m_dot*cp*(T3-T2); 
Energy_q_rate_loss_heatingair = 0;
%% Entropy Production in Heating Air System
s3_s2 = cp * log(T3 / T2) - R * log(P3 / P2); % Entropy change (J/kg·K)
s3_s0 = cp * log(T3 / T0) - R * log(P3 / P0); % Entropy change (J/kg·K)

sigma_dot_p_heatingair= -(Energy_q_rate_in_heatingair/T_flame) + m_dot * s3_s2; % TH_hp_actual% Entropy production rate (W/K)
Exg_d_rate_heatingair= T0 * sigma_dot_p_heatingair; % Exergy destroyed (W)

%% Exergy Analysis for Heating Air SYstem

% Exergy flow
exg_f_in_heatingair= exg_f_out_fan; % Specific exergy at Heating by Furnace heating inlet (J/kg)
exg_f_out_heatingair= cp * (T3 - T0) - T0 * s3_s0; % Specific exergy at outlet (J/kg)
Exg_f_rate_in_heatingair= m_dot * exg_f_in_heatingair; % Exergy inflow (W)
Exg_f_rate_out_heatingair= m_dot * exg_f_out_heatingair; % Exergy outflow (W)

% Exergy in/out
Exg_q_rate_in_heatingair= (1-T0/T_flame)*Q_dot_heatingair; %Q_dot_actual_heating;
Exg_w_rate_in_heatingair= 0;
Exg_rate_in_heatingair=Exg_q_rate_in_heatingair+Exg_w_rate_in_heatingair;

Exg_q_rate_out_heatingair = 0;
Exg_w_rate_out_heatingair=0;
Exg_rate_out_heatingair=Exg_q_rate_out_heatingair+Exg_w_rate_out_heatingair;
Exg_q_rate_loss_heatingair = 0;

Exg_d_rate_heatingair_exgbalance= (Exg_q_rate_in_heatingair-Exg_q_rate_out_heatingair) +(Exg_w_rate_in_heatingair-Exg_w_rate_out_heatingair)...
    +(-Exg_q_rate_loss_heatingair)+(Exg_f_rate_in_heatingair -Exg_f_rate_out_heatingair);

Exg_utilized_rate_heatingair =Exg_f_rate_out_heatingair-Exg_f_rate_in_heatingair;
%  Efficiency

eff_energy_heatingair= Energy_utilized_rate_heatingair/Energy_rate_in_heatingair;
eff_exergy_heatingair= Exg_utilized_rate_heatingair/Exg_rate_in_heatingair;

%% Print Results
fprintf('\nExergy Destroyed Rate in Heating Air Sub-system (W): %.2f\n', Exg_d_rate_heatingair);
fprintf('Exergy Destroyed Rate in Heating Air Sub-system (checked by Exergy Balnace) (W): %.2f\n',Exg_d_rate_heatingair_exgbalance)

%fprintf('Energy Efficiencey of Heating Air Sub-system(%%): %.2f\n', eff_energy_heatingair*100);
%fprintf('Exergy Efficiencey of Heating AIr Sub-system (%%): %.2f\n', eff_exergy_heatingair*100);
end