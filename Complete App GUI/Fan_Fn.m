function [T2, W_dot_ideal_fan,...
    Energy_w_rate_in_fan, Energy_q_rate_in_fan, Energy_rate_in_fan, ...
    Energy_w_rate_out_fan, Energy_q_rate_out_fan, Energy_rate_out_fan, Energy_q_rate_loss_fan,...
    ...
    Exg_w_rate_in_fan, Exg_q_rate_in_fan, Exg_rate_in_fan, ...
    Exg_w_rate_out_fan, Exg_q_rate_out_fan, Exg_rate_out_fan, Exg_q_rate_loss_fan, ...
    ...
    exg_f_in_fan, exg_f_out_fan, Exg_f_rate_in_fan, Exg_f_rate_out_fan, ...
    sigma_dot_p_fan, Exg_d_rate_fan, Exg_d_rate_fan_exgbalance, ...
    ...
    Energy_utilized_rate_fan, Exg_utilized_rate_fan, eff_exergy_fan, eff_energy_fan] = Fan_Fn(m_dot, cp, R, gamma, T0, T1, P0, P1, P2, P_design,V_dot,eta_fan)
%% Premlinminary Analysis - Fan Power and States
delta_P = P2 - P1; % Pressure difference across fan (Pa)
W_dot_ideal_fan = V_dot * delta_P; % Ideal fan power (W)

% Actual Fan Power
W_dot_actual_fan = W_dot_ideal_fan / eta_fan; % Actual fan power (W)

% Fan Outlet Temperature
T2_ideal = T1 * (P_design / P1)^((gamma - 1) / gamma); % Ideal outlet temperature (K)
T2_actual=T1+(W_dot_actual_fan/m_dot/cp); % Actual outlet temperature (K)
T2=T2_actual;

%% Energy Analysis for Fan
Energy_q_rate_in_fan = 0;
Energy_w_rate_in_fan = W_dot_actual_fan; % Fan energy input (W)
Energy_rate_in_fan = Energy_q_rate_in_fan+Energy_w_rate_in_fan;

Energy_q_rate_out_fan = 0;
Energy_w_rate_out_fan = 0;
Energy_rate_out_fan = Energy_w_rate_out_fan + Energy_q_rate_out_fan;
Energy_q_rate_loss_fan = 0;

Energy_utilized_rate_fan = m_dot*cp*(T2-T1);

%% Entropy Production in Fan

s1_s0 = cp * log(T1 / T0) - R * log(P1 / P0); % Only for inlet of first component, for next ones, it will use the precedent component

s2_s1 = cp * log(T2 / T1) - R * log(P2 / P1); % Entropy change (J/kg·K)
s2_s0 = cp * log(T2 / T0) - R * log(P2 / P0); % Entropy change (J/kg·K)

sigma_dot_p_fan = m_dot * s2_s1; % No - Q in and Q out; Entropy production rate (W/K)
Exg_d_rate_fan = T0 * sigma_dot_p_fan; % Exergy destroyed (W)

%% Exergy Analysis for Fan

% Exergy flow
exg_f_in_fan = cp * (T1 - T0) - T0 * s1_s0; % Inlet specific exergy (J/kg)
exg_f_out_fan = cp * (T2 - T0) - T0 * s2_s0; % Outlet specific exergy (J/kg)
Exg_f_rate_in_fan = m_dot * exg_f_in_fan; % Exergy inflow (W)
Exg_f_rate_out_fan = m_dot * exg_f_out_fan; % Exergy outflow (W)

% Exergy in/out
Exg_q_rate_in_fan = 0;
Exg_w_rate_in_fan = W_dot_actual_fan;
Exg_rate_in_fan = Exg_q_rate_in_fan+ Exg_w_rate_in_fan; % Fan work exergy input (W)

Exg_q_rate_out_fan = 0; % NOT loss
Exg_w_rate_out_fan = 0; 
Exg_rate_out_fan = Exg_q_rate_out_fan+ Exg_w_rate_out_fan; % Fan work exergy input (W)
Exg_q_rate_loss_fan = 0;

Exg_d_rate_fan_exgbalance = (Exg_q_rate_in_fan-Exg_q_rate_out_fan) +(Exg_w_rate_in_fan-Exg_w_rate_out_fan)+...
    +(-Exg_q_rate_loss_fan)+(Exg_f_rate_in_fan -Exg_f_rate_out_fan);
Exg_utilized_rate_fan = Exg_f_rate_out_fan-Exg_f_rate_in_fan;

%% Energy/ Exergy Efficiencey

eff_exergy_fan = Exg_utilized_rate_fan/Exg_rate_in_fan;
eff_energy_fan = Energy_utilized_rate_fan/Energy_rate_in_fan; %eff_energy_fan = W_dot_ideal_fan/W_dot_actual_fan;

%% Display Results
% Energy/Exergy Results -Fan
disp('Fan Results:');
fprintf('Ideal Fan Power (W): %.2f\n', W_dot_ideal_fan);
fprintf('Actual Fan Power (W): %.2f\n', W_dot_actual_fan);

% 
fprintf('Exergy Destroyed rate in Fan (W): %.2f\n', Exg_d_rate_fan);
fprintf('Exergy Destroyed rate in Fan - Exergy balance (W): %.2f\n', Exg_d_rate_fan_exgbalance );

%fprintf('Energy Efficiencey of Fan (%%): %.2f\n', eff_energy_fan*100);
%fprintf('Exergy Efficiencey of Fan (%%): %.2f\n', eff_exergy_fan*100);

end