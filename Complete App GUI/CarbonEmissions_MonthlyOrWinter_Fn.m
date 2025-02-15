function [electricity_geneartaion_kWh_period_analyzed,electricity_fossil_CO2_emissions_period_analyzed,...
    natural_gas_CO2_emissions_period_analyzed, total_CO2_emissions_period_analyzed,Carbon_Emissions_Table_Period]=...
    CarbonEmissions_MonthlyOrWinter_Fn(name_period_analyzed,Mol_CO2,Mol_CH4,HHV_ng,eff_power_plant_ng,...
    density_ng,renewable_share_electricity,transmission_loss_electricity, ...
    Energy_in_fan_kWh_period_analyzed, natural_gas_volume_m3_period_analyzed)
%% CO2 Production % CO2 Emission Analysis

m_ratio_CO2_CH4 = Mol_CO2/Mol_CH4; %mol_CO2 / Mol_CH4 ==> Mol_CO2/ Mol_Ng
m_ratio_CO2_Ng = 1 * m_ratio_CO2_CH4; % same ratio

emission_factor_fossil_electricity = m_ratio_CO2_Ng/(HHV_ng*eff_power_plant_ng)*3.6; 
% kg CO2/kWh (fossil fuel contribution for electricity generation in Boston)
emission_factor_burning_ng = m_ratio_CO2_Ng * density_ng ; %[kgCO2/m3Ng] for natural gas combustion
% 2.74 kgCO2/kgNg * ~kgNg/m3Ng =>  [kgCO2/m3Ng] for natural gas combustion

fossil_share_electricity = 1 - renewable_share_electricity; % Fossil fuel contribution

% Effective electricity usage (including transmission losses)
electricity_geneartaion_kWh_period_analyzed = Energy_in_fan_kWh_period_analyzed / (1 - transmission_loss_electricity);

%% CO2 emission
% CO2 emissions from fossil fuel-based electricity
electricity_fossil_CO2_emissions_period_analyzed = electricity_geneartaion_kWh_period_analyzed * fossil_share_electricity * emission_factor_fossil_electricity;

% Total CO2 emissions from natural gas burning in furanace
natural_gas_CO2_emissions_period_analyzed = natural_gas_volume_m3_period_analyzed * emission_factor_burning_ng;

total_CO2_emissions_period_analyzed = electricity_fossil_CO2_emissions_period_analyzed + natural_gas_CO2_emissions_period_analyzed;

%% Display CO2 Emission Results
%{
fprintf('CO2 Emission Analysis Results in "%s"\n',name_period_analyzed);
fprintf('Electricity CO2 Emissions including Plant Production (kg): %.2f\n', electricity_fossil_CO2_emissions_period_analyzed);
fprintf('Natural Gas CO2 Emissions for Furnace Heating by burning  (kg): %.2f\n', natural_gas_CO2_emissions_period_analyzed);
fprintf('Total CO2 Emissions (kg): %.2f\n', total_CO2_emissions_period_analyzed);
%}
%% Updated Cost and Emission Summary Table
Component = {'Fan'; 'Furnace'};
Energy_Input_kWh = [Energy_in_fan_kWh_period_analyzed; 0];
Natural_Gas_Input_m3 = [0; natural_gas_volume_m3_period_analyzed];
CO2_Emissions_kg = [electricity_fossil_CO2_emissions_period_analyzed; natural_gas_CO2_emissions_period_analyzed];

Carbon_Emissions_Table_Period = table(Component, Energy_Input_kWh, Natural_Gas_Input_m3, CO2_Emissions_kg, ...
    VariableNames={'Source System','Energy Input Units [kWh]','Fuel Input Volume [m3]','C02 emisssions [kg]'});
fprintf('\nTable 7: Energy Source and CO2 Emission Summary in %s\n',name_period_analyzed);
disp(Carbon_Emissions_Table_Period);
end