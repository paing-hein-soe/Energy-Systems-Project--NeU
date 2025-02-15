function [State_Table, states, Specific_Energy, Carried_Energy_Rate, Specific_Flow_Exergy, Flow_Exergy_Rate] = States_Fn(m_dot,cp,R,T0, T1, T2, T3, T4, P0, P1, P2, P3, P4)

specific_flow_exergy = @(T, P) cp * (T - T0) - T0 * (cp * log(T / T0) - R * log(P / P0));

%% Define States Dynamically
states = {
    struct('label', '1', 'T', T1, 'P', P1), ...
    struct('label', '2', 'T', T2, 'P', P2), ...
    struct('label', '3', 'T', T3, 'P', P3), ...
    struct('label', '4', 'T', T4, 'P', P4)
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
State_Label = {states{1}.label; states{2}.label; states{3}.label; states{4}.label};
Description = {'Fan Inlet'; 'Fan Outlet'; 'Heating by Furnace- Outlet'; 'Hospital Discharge'};
Temperature_K = [states{1}.T; states{2}.T; states{3}.T; states{4}.T];
Pressure_Pa = [states{1}.P; states{2}.P; states{3}.P; states{4}.P];

State_Table = table(State_Label, Description, Temperature_K, Pressure_Pa, ...
    Carried_Energy_Rate, Specific_Energy, Flow_Exergy_Rate, Specific_Flow_Exergy, ...
    'VariableNames', {'State_Label', 'Description', 'Temperature [K]', 'Pressure [Pa]', ...
    'Carried_Energy_Rate [W]', 'Specific_Energy [J/kg]', 'Flow_Exergy_Rate [W]', 'Specific_Flow_Exergy [J/kg]'});

disp('Table 1: State Properties and Energy/Exergy Analysis');
disp(State_Table);
end
