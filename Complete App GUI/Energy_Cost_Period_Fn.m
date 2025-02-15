function [specific_volume_ng, natural_gas_energy_MJ_per_m3, Energy_in_fan_MJ_period, Energy_in_fan_kWh_period, electricity_cost_fan_period,...
    Energy_in_heatingwhole_MJ_period, natural_gas_volume_m3_period, fuel_cost_furnace_period]=Energy_Cost_Period_Fn(density_ng,...
    HHV_ng, Energy_in_fan_period,Energy_in_heatingwhole_period,electricity_cost_per_kWh,natural_gas_cost_per_m3)
    %% Cost Analysis -Electricity & Fuel Furance
    specific_volume_ng = 1/density_ng ; %m^3/kg
    natural_gas_energy_MJ_per_m3 = HHV_ng/specific_volume_ng; % MJ/m³
    
    % Step 1: Electricity Cost for Fan
    Energy_in_fan_MJ_period = Energy_in_fan_period/1e6; % MJ/winterseason (from Table 4)
    Energy_in_fan_kWh_period = Energy_in_fan_MJ_period / 3.6; % Convert MJ to kWh
    electricity_cost_fan_period = Energy_in_fan_kWh_period * electricity_cost_per_kWh; % Total electricity cost for fan ($)
    
    % Step 2: Fuel Cost for Heating by Furnace
    Energy_in_heatingwhole_MJ_period = Energy_in_heatingwhole_period/1e6; % MJ (from Table 4, heat input for furnace)
    natural_gas_volume_m3_period = Energy_in_heatingwhole_MJ_period / natural_gas_energy_MJ_per_m3; % Volume of natural gas required (m³)
    fuel_cost_furnace_period = natural_gas_volume_m3_period * natural_gas_cost_per_m3; % Total fuel cost for furnace ($)

end