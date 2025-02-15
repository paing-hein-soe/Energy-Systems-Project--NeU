function [State_Table, states, Specific_Energy, Carried_Energy_Rate, Specific_Flow_Exergy, Flow_Exergy_Rate] = States_Option4_Fn(m_dot,cp,R,T0, T1, T2, T2A, T2B, T3, T4, T4A, P0, P1, P2, P2A, P2B, P3, P4, P4A)

specific_flow_exergy = @(T, P) cp * (T - T0) - T0 * (cp * log(T / T0) - R * log(P / P0));

%% Define States Dynamically
states = {
    struct('label', '1', 'T', T1, 'P', P1), ...
    struct('label', '2', 'T', T2, 'P', P2), ...
    struct('label', '2A', 'T', T2A, 'P', P2A), ...
    struct('label', '2B', 'T', T2B, 'P', P2B), ...
    struct('label', '3', 'T', T3, 'P', P3), ...
    struct('label', '4', 'T', T4, 'P', P4)...
    struct('label', '4A', 'T', T4A, 'P', P4A)
};

% Initialize Arrays for Storing Results
Specific_Energy = zeros(length(states), 1);
Carried_Energy_Rate = zeros(length(states), 1);
Specific_Flow_Exergy = zeros(length(states), 1);
Flow_Exergy_Rate = zeros(length(states), 1);

% Calculate Values for Each State
for i = 1:length(states)
    T = states{i}.T;
    P = states{i}.P;
    h = cp * T; % Specific energy (J/kg)
    carried_energy_rate = m_dot * h; % Carried energy rate (W)
    e_f = specific_flow_exergy(T, P); % Specific flow exergy (J/kg)
    flow_exergy_rate = m_dot * e_f; % Flow exergy rate (W)
    
    % Store Results
    Specific_Energy(i) = h;
    Carried_Energy_Rate(i) = carried_energy_rate;
    Specific_Flow_Exergy(i) = e_f;
    Flow_Exergy_Rate(i) = flow_exergy_rate;
end

%% Generate Table 1 - State Properties
State_Label = {states{1}.label; states{2}.label; states{3}.label; states{4}.label; states{5}.label;states{6}.label;states{7}.label};
Description = {'Fan Inlet'; 'Fan Outlet';'2A (After Heat Transfer from RegHE )';'2B (After heating by HE)'; '3 (After heating by HP)'; '4-Hospital Discharge';'4A-Outlet of RegHE'};
Temperature_K = [states{1}.T; states{2}.T; states{3}.T; states{4}.T;states{5}.T;states{6}.T;states{7}.T];
Pressure_Pa = [states{1}.P; states{2}.P; states{3}.P; states{4}.P; states{5}.P;states{6}.P;states{7}.P];

State_Table = table(State_Label, Description, Temperature_K, Pressure_Pa, ...
    Carried_Energy_Rate, Specific_Energy, Flow_Exergy_Rate, Specific_Flow_Exergy, ...
    'VariableNames', {'State_Label', 'Description', 'Temperature [K]', 'Pressure [Pa]', ...
    'Carried_Energy_Rate [W]', 'Specific_Energy [J/kg]', 'Flow_Exergy_Rate [W]', 'Specific_Flow_Exergy [J/kg]'});

disp('Table 1: State Properties and Energy/Exergy Analysis');
disp(State_Table);
end
