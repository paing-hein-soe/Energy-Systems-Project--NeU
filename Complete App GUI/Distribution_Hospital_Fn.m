function  [Energy_w_rate_in_hospital, Energy_q_rate_in_hospital, Energy_rate_in_hospital, ...
    Energy_w_rate_out_hospital, Energy_q_rate_out_hospital, Energy_rate_out_hospital, Energy_q_rate_loss_hospital,...
    ...
    Exg_w_rate_in_hospital, Exg_q_rate_in_hospital, Exg_rate_in_hospital, ...
    Exg_w_rate_out_hospital, Exg_q_rate_out_hospital, Exg_rate_out_hospital, Exg_q_rate_loss_hospital, ...
    ...
    exg_f_in_hospital, exg_f_out_hospital, Exg_f_rate_in_hospital, Exg_f_rate_out_hospital, ...
    sigma_dot_p_hospital, Exg_d_rate_hospital, Exg_d_rate_hospital_exgbalance, ...
    ...
    Energy_utilized_rate_hospital, Exg_utilized_rate_hospital, eff_exergy_hospital, eff_energy_hospital]=Distribution_Hospital_Fn(m_dot, cp, R, T0,T3,T4, P0, P3, P4)

%% Premlinminary Analysis
Q_dot_out_hospital= m_dot * cp * (T3 - T4); % Heat transfer (W)
T_boundary_hospital = T3; % Boundary temperature for heat loss (K)

%% Energy Analysis for Hospital Distribution
Energy_q_rate_in_hospital= 0;
Energy_w_rate_in_hospital= 0;
Energy_rate_in_hospital = m_dot * cp * (T3 - T4)+Energy_q_rate_in_hospital+ Energy_w_rate_in_hospital;

Energy_q_rate_out_hospital= Q_dot_out_hospital;
Energy_w_rate_out_hospital=0;
Energy_rate_out_hospital= Energy_q_rate_out_hospital+Energy_w_rate_out_hospital;
Energy_q_rate_loss_hospital = 0; % NO losses assumptions

Energy_utilized_rate_hospital = Q_dot_out_hospital;

%% Entropy Production in Hospital
s3_s0 =cp * log(T3 / T0) - R * log(P3 / P0);
s4_s3 = cp * log(T4 / T3) - R * log(P4 / P3); % Entropy change (J/kg·K)
s4_s0 = cp * log(T4 / T0) - R * log(P4 / P0);

sigma_dot_p_hospital = -(-Q_dot_out_hospital / T_boundary_hospital) + (m_dot * s4_s3); % Entropy production rate (W/K)
Exg_d_rate_hospital = T0 * sigma_dot_p_hospital;

%% Exergy Analysis in Hospital Air Distribution
exg_f_in_hospital = cp * (T3 - T0) - T0 * s3_s0 ; %exg_f_out_heatingwhole;
exg_f_out_hospital = cp * (T4 - T0) - T0 * s4_s0 ; % Specific exergy (J/kg)

Exg_f_rate_in_hospital = m_dot * exg_f_in_hospital; %Exg_f_rate_out_heatingwhole; 
Exg_f_rate_out_hospital = m_dot * exg_f_out_hospital; % Exergy outflow (W)

% Exergy in/out
Exg_q_rate_in_hospital = 0 ;
Exg_w_rate_in_hospital = 0 ;
Exg_rate_in_hospital = (Exg_f_rate_in_hospital-Exg_f_rate_out_hospital); %Take flow rate as input

Exg_q_rate_out_hospital = (1 - T0 / T_boundary_hospital) * (Q_dot_out_hospital); % Exergy by heat transfer
Exg_w_rate_out_hospital = 0;
Exg_rate_out_hospital = Exg_q_rate_out_hospital+Exg_w_rate_out_hospital;
Exg_q_rate_loss_hospital = 0;

Exg_d_rate_hospital_exgbalance = (Exg_q_rate_in_hospital-Exg_q_rate_out_hospital) +(Exg_w_rate_in_hospital-Exg_w_rate_out_hospital)+...
    +(-Exg_q_rate_loss_hospital)+(Exg_f_rate_in_hospital -Exg_f_rate_out_hospital);

Exg_utilized_rate_hospital =  Exg_f_rate_in_hospital -Exg_f_rate_out_hospital; 

%% Energy/ Exergy Efficiencey
eff_energy_hospital= Energy_utilized_rate_hospital/ Energy_rate_in_hospital;
eff_exergy_hospital= Exg_utilized_rate_hospital/Exg_rate_in_hospital;

%% Display results
fprintf('\nHospital Air Distribution Results:');
fprintf('Heat Transfer in Hospital Path (W): %.2f\n', Q_dot_out_hospital);

fprintf('Exergy Destroyed Rate in Hospital Path (W): %.2f\n', Exg_d_rate_hospital);
fprintf('Exergy Destroyed Rate in Hospital Path (checked by exergy balance) (W): %.2f\n', Exg_d_rate_hospital_exgbalance);

%fprintf('Hospital Air Distribution - First-Law (Energy) Efficiency: %.2f%%\n', eff_energy_hospital * 100);
%fprintf('Hospital Air Distribution - Second-Law (Exergy) Efficiency: %.2f%%\n', eff_exergy_hospital * 100);

end