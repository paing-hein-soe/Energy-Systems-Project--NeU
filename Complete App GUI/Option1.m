function [Cost_Table, Carbon_Emissions_Table_winter, plotData]=Option1(T0_list_input,P0_input,T_design_input,P_design_input,T_discharge_input,P_disccharge_input)
% MATLAB Code for Fan and Heyating by Furnace Calculations

clc;  close all;
%% Temperature Change Effect

T0_list = T0_list_input;%-10:2:20; %8;% T0 = 0; % Ambient temperature (°C)
P0_list = P0_input *ones(size(T0_list)); %P1 = 1; % Design pressure (bar)

T_design_list = T_design_input* ones(size(T0_list)); % Base: T_design = 36; % Design temperature (°C)
P_design_list = P_design_input*ones(size(T0_list)); %P_design = 1.1; % Design pressure (bar)

T4_list = T_discharge_input*ones(size(T0_list));%24* ones(size(T0_list)); % Base: T_design = 36; % Design temperature (°C)
P4_list = P_disccharge_input*ones(size(T0_list)); %*ones(size(T0_list)); %P_design = 1.1; % Design pressure (bar)

%% Initializations For Temeperature Change Effect 
    Energy_in_fan_kWh_winter_list_Teffect=zeros(length(T0_list),1);
    electricity_cost_fan_winter_list_Teffect = zeros(length(T0_list),1);

    natural_gas_volume_m3_winter_list_Teffect =zeros(length(T0_list),1);
    fuel_cost_furnace_winter_list_Teffect = zeros(length(T0_list),1);

    total_cost_withCapital_period_winter_list_Teffect = zeros(length(T0_list),1);
    
    % CO2 Emissions
    electricity_fossil_CO2_emissions_winter_list_Teffect= zeros(length(T0_list),1);
    natural_gas_CO2_emissions_winter_list_Teffect = zeros(length(T0_list),1);
    
    total_CO2_emissions_winter_list_Teffect = zeros(length(T0_list),1);

    % Exergy destruction and Exe
    Energetic_efficiency_winter_list_Teffect = zeros(length(T0_list),3); %for fan, heating and furnace
    Exergetic_efficiency_winter_list_Teffect =  zeros(length(T0_list),3); %for fan, heating and furnace

for tt= 1:length(T0_list)
        
    %% Part I - Design Inputs - 
    % Ambient Conditions
    T0 = T0_list(tt); % Ambient temperature (°C)
    P0 = P0_list(tt); % Ambient pressure (bar)
    
    % Design Conditions
    T_design = T_design_list(tt); % Design temperature (°C)
    P_design = P_design_list(tt); % Design pressure (bar)
    
    % Discharge State
    T4 = T4_list(tt); % Discharge temperature (°C)
    P4 = P4_list(tt); % Discharge pressure (bar)
    % Needs - Unit Conversion for T [K] and P [Pa]
    
    %% Defined States
    T1 = T0; T3 = T_design;
    P1 = P0; P2 = P_design; P3 = P_design;  
    T_flame = 1956.85; % Flame Temperaute of furnace in (°C)

    %% Unit Converstions
    T0 = T0 + 273.15; T_design = T_design + 273.15; T4 = T4 + 273.15; T_flame = T_flame + 273.15; % C to Kelvin
    T1 = T1 + 273.15; T3 = T3 + 273.15; 
    P0 = P0*1e5; P_design = P_design * 1e5; P4 = P4 * 1e5; % bar to Pa
    P1 = P1*1e5; P2 = P2*1e5; P3 = P3*1e5; 
    %% Part (II) Seasonal Energy/Exergy Analysis
    months_design = {'November', 'December', 'January', 'February', 'March', 'April'};% Winter months
    [days_per_month, hours_per_month, seconds_per_month, total_days, total_hours, total_seconds]=duration_fn(months_design);
    %dt_analysis = total_seconds; % <== Can change eg. seconds per month ==> It can be looped on T changes over dt period
    analysis_period = 'Winter';
    
    %% Part II-Design Inputs -Seasonal Energy/Exergy Analysis 
    %% Default Parameters - Except Ambient and Design Temperatures
    % Constants for CO2 Emission Factors
    Mol_CH4 = 16.04; % g/mol; 
    Mol_CO2 = 44.01; % g/mol; 
    
    % Operation Conditions
    V_building = 10000; % Building volume (m^3)
    N_changes = 12; % Air changes per hour
    
    % Fan Performance
    eta_fan = 0.6; % Fan efficiency
    
    % Furnace Performance
    eta_furnace= 0.9; % Heating by Furnace efficiency
    
    %% Electricity & Fuel 
    electricity_cost_per_kWh = 0.2102; % $/kWh * for commercial building
    natural_gas_cost_per_m3 = 0.82177; % $/m³
    
    HHV_ng = 55.5; %MJ/kg
    density_ng = 0.668; %kg/m3
    
    eff_power_plant_ng = 0.5; % 50%
    renewable_share_electricity = 0.32; %  32 percent renewable and remaining fossil fuels (Ng)
    transmission_loss_electricity = 0.07; % 7% loss during electricity transmission
    
    %% Fan/ Furnace Finance
    % Capital Costs (Assumed Values, Will update with More Research)
    capital_cost_fan = 10000; % $ (Assume purchase and installation cost)
    capital_cost_furnace = 15000; % $ (Assume purchase and installation cost)
    
    %
    capital_cost_he = 0; % $100,000 
    captial_cost_hp = 0; % $50,000
    capital_cost_regHE = 0; %$10,000

    % Financing Period and Annual Interest Rate
    financing_period_years = 10; % Financing over 10 years
    annual_interest_rate = 0.12; % 12% annual interest rate
    
    operator_full_cost_per_year = 95000; % $/year
    operator_allocation = 0.2; % 20% of operator time
    
    %% Gas Constants
    cp = 1005; % Specific heat of air at constant pressure (J/kg·K)
    R = 287;   % Specific gas constant for air (J/kg·K)
    gamma = 1.4; % Ratio of specific heats (cp/cv)
    
    %% Checking Monthly energy loops
    switch analysis_period
        case 'Hourly'
            dt_analysis_period = 60*ones(size(seconds_per_month));
            monthly_loops = hours_per_month;
        case 'Bi-hourly'
            dt_analysis_period = 120*ones(size(seconds_per_month));
            monthly_loops = hours_per_month/2; 
        case 'Monthly'
            dt_analysis_period = seconds_per_month;
            monthly_loops = ones(size(months_design));
        case 'Winter'
            dt_analysis_period = total_seconds;
            monthly_loops = 1;
    end
    %% Initializations for Monthly Data 
    if ~strcmpi(analysis_period,'Winter')
        Energy_in_fan_MJ_monthly = zeros(length(months_design),1);
        Energy_in_fan_kWh_monthly = zeros(length(months_design),1);
        electricity_cost_fan_monthly = zeros(length(months_design),1);
        
        Energy_in_heatingwhole_MJ_monthly = zeros(length(months_design),1);
        natural_gas_volume_m3_monthly = zeros(length(months_design),1);
        fuel_cost_furnace_monthly = zeros(length(months_design),1);
    end
    %% 
    for ii = 1:length(monthly_loops)
        Energy_in_fan_MJ_period_list = zeros(length(monthly_loops(ii)),1); 
        Energy_in_fan_kWh_period_list = zeros(length(monthly_loops(ii)),1);
        electricity_cost_fan_period_list = zeros(length(monthly_loops(ii)),1);
        
        Energy_in_heatingwhole_MJ_period_list =  zeros(length(monthly_loops(ii)),1);
        natural_gas_volume_m3_period_list =  zeros(length(monthly_loops(ii)),1);
        fuel_cost_furnace_period_list=  zeros(length(monthly_loops(ii)),1);
    
        dt_analysis = dt_analysis_period(ii);
        for jj = 1: monthly_loops(ii)
    

    %% Analysis Start 
    %% Fan Analysis
    % Initial calculations 
    V_dot = (N_changes * V_building) / 3600; % m^3/s
    v_outside = R * T1 / P1; % Specific volume at outside conditions (m^3/kg)
    m_dot = V_dot / v_outside; % Mass flow rate (kg/s)
    
    [T2, W_dot_ideal_fan,...
        Energy_w_rate_in_fan, Energy_q_rate_in_fan, Energy_rate_in_fan, ...
        Energy_w_rate_out_fan, Energy_q_rate_out_fan, Energy_rate_out_fan, Energy_q_rate_loss_fan,...
        ...
        Exg_w_rate_in_fan, Exg_q_rate_in_fan, Exg_rate_in_fan, ...
        Exg_w_rate_out_fan, Exg_q_rate_out_fan, Exg_rate_out_fan, Exg_q_rate_loss_fan, ...
        ...
        exg_f_in_fan, exg_f_out_fan, Exg_f_rate_in_fan, Exg_f_rate_out_fan, ...
        sigma_dot_p_fan, Exg_d_rate_fan, Exg_d_rate_fan_exgbalance, ...
        ...
        Energy_utilized_rate_fan, Exg_utilized_rate_fan, eff_exergy_fan, eff_energy_fan] = Fan_Fn(m_dot, cp, R, gamma, T0, T1, P0, P1, P2, P_design,V_dot,eta_fan);
    
    %% Heating Air - Analysis
    Q_dot_heatingair= m_dot * cp * (T_design - T2); % Heat added to air (W)% Heat Required for Heating by Furnace
    
    % --> Function For Heating Air 
    [Energy_w_rate_in_heatingair, Energy_q_rate_in_heatingair, Energy_rate_in_heatingair, ...
        Energy_w_rate_out_heatingair, Energy_q_rate_out_heatingair, Energy_rate_out_heatingair, Energy_q_rate_loss_heatingair,...
        ...
        Exg_w_rate_in_heatingair, Exg_q_rate_in_heatingair, Exg_rate_in_heatingair, ...
        Exg_w_rate_out_heatingair, Exg_q_rate_out_heatingair, Exg_rate_out_heatingair, Exg_q_rate_loss_heatingair, ...
        ...
        exg_f_in_heatingair, exg_f_out_heatingair, Exg_f_rate_in_heatingair, Exg_f_rate_out_heatingair, ...
        sigma_dot_p_heatingair, Exg_d_rate_heatingair, Exg_d_rate_heatingair_exgbalance, ...
        ...
        Energy_utilized_rate_heatingair, Exg_utilized_rate_heatingair, eff_exergy_heatingair, eff_energy_heatingair]=HeatingAir_Fn...
        (m_dot, cp, R, T0, T2, T3, P0, P2, P3, Q_dot_heatingair, T_flame, exg_f_out_fan );
    
    
    %% Option 1 - Heating by Furnace Analysis :Direct Heating by Furnace
    [Energy_w_rate_in_furnace, Energy_q_rate_in_furnace, Energy_rate_in_furnace, ...
        Energy_w_rate_out_furnace, Energy_q_rate_out_furnace, Energy_rate_out_furnace, Energy_q_rate_loss_furnace,...
        ...
        Exg_w_rate_in_furnace, Exg_q_rate_in_furnace, Exg_rate_in_furnace, ...
        Exg_w_rate_out_furnace, Exg_q_rate_out_furnace, Exg_rate_out_furnace, Exg_q_rate_loss_furnace, ...
        ...
        ...
        sigma_dot_p_furnace, Exg_d_rate_furnace, Exg_d_rate_furnace_exgbalance, ...
        ...
        Energy_utilized_rate_furnace, Exg_utilized_rate_furnace, eff_exergy_furnace, eff_energy_furnace]=Furnace_Direct_Option1_Fn...
    (T0, T_flame, eta_furnace, Energy_q_rate_in_heatingair, Energy_w_rate_in_heatingair, Exg_q_rate_in_heatingair, Exg_w_rate_in_heatingair);
    
    %% Option 1- Whole heating analysis
    
    [Energy_w_rate_in_heatingwhole, Energy_q_rate_in_heatingwhole, Energy_rate_in_heatingwhole,...
        Energy_w_rate_out_heatingwhole, Energy_q_rate_out_heatingwhole, Energy_rate_out_heatingwhole, Energy_q_rate_loss_heatingwhole,...
        ...
        Exg_w_rate_in_heatingwhole, Exg_q_rate_in_heatingwhole, Exg_rate_in_heatingwhole, ...
        Exg_w_rate_out_heatingwhole, Exg_q_rate_out_heatingwhole, Exg_rate_out_heatingwhole, Exg_q_rate_loss_heatingwhole, ...
        ...
        exg_f_in_heatingwhole, exg_f_out_heatingwhole, Exg_f_rate_in_heatingwhole, Exg_f_rate_out_heatingwhole, ...
        sigma_dot_p_heatingwhole, Exg_d_rate_heatingwhole, Exg_d_rate_heatingwhole_exgbalance, ...
        ...
        Energy_utilized_rate_heatingwhole, Exg_utilized_rate_heatingwhole, eff_exergy_heatingwhole, eff_energy_heatingwhole]=HeatingWhole_Direct_Option1_Fn...
        ( Energy_q_rate_in_furnace, Energy_w_rate_in_furnace,Energy_q_rate_out_heatingair,Energy_w_rate_out_heatingair, ...
        Energy_utilized_rate_heatingair,Energy_q_rate_loss_heatingair,Energy_q_rate_loss_furnace,...
        ...
        sigma_dot_p_furnace,sigma_dot_p_heatingair,Exg_d_rate_furnace,Exg_d_rate_heatingair,...
        ...
        exg_f_in_heatingair,exg_f_out_heatingair,Exg_f_rate_in_heatingair,Exg_f_rate_out_heatingair,...
        Exg_q_rate_in_furnace, Exg_w_rate_in_furnace,Exg_q_rate_out_heatingair, Exg_w_rate_out_heatingair,...
        Exg_utilized_rate_heatingair, Exg_q_rate_loss_heatingair, Exg_q_rate_loss_furnace);
    
    %% Hospital Air Distribution Air Analysis
    %Q_dot_out_hospital= m_dot * cp * (T3 - T4); % Heat transfer (W)
    [Energy_w_rate_in_hospital, Energy_q_rate_in_hospital, Energy_rate_in_hospital, ...
        Energy_w_rate_out_hospital, Energy_q_rate_out_hospital, Energy_rate_out_hospital, Energy_q_rate_loss_hospital,...
        ...
        Exg_w_rate_in_hospital, Exg_q_rate_in_hospital, Exg_rate_in_hospital, ...
        Exg_w_rate_out_hospital, Exg_q_rate_out_hospital, Exg_rate_out_hospital, Exg_q_rate_loss_hospital, ...
        ...
        exg_f_in_hospital, exg_f_out_hospital, Exg_f_rate_in_hospital, Exg_f_rate_out_hospital, ...
        sigma_dot_p_hospital, Exg_d_rate_hospital, Exg_d_rate_hospital_exgbalance, ...
        ...
        Energy_utilized_rate_hospital, Exg_utilized_rate_hospital, eff_exergy_hospital, eff_energy_hospital]=Distribution_Hospital_Fn(m_dot, cp, R, T0,T3,T4, P0, P3, P4);
    
    %% Overall System - Energy/Exergy/Efficiency Analysis 
    
    [Energy_rate_in_system,Exg_rate_disposed_system, Energy_utilized_rate_system,...
        Exg_rate_in_system, Exg_utilized_rate_system, Exg_d_rate_system, Exg_q_rate_loss_system,...
        eff_energy_system, eff_exergy_system]= Overall_EneExeEff_Fn(...
        Energy_rate_in_fan, Energy_rate_in_heatingwhole,Energy_utilized_rate_hospital,...
        Exg_rate_in_fan, Exg_rate_in_heatingwhole,Exg_utilized_rate_hospital,...
        Exg_d_rate_fan, Exg_d_rate_heatingwhole, Exg_d_rate_hospital,...
        Exg_q_rate_loss_fan, Exg_q_rate_loss_heatingwhole, Exg_q_rate_loss_hospital,Exg_f_rate_out_hospital);
    
    %% Overall system - Exergy Balance Check Table
    [Exg_rate_in_system_check,Exg_rate_disposed_system_check,exergy_sources_table,exergy_disposition_table]=...
        Exg_Balance_Table_Option1_Fn(m_dot,Exg_rate_in_system,Exg_f_rate_in_fan,Exg_w_rate_in_fan,Exg_q_rate_in_heatingwhole,...
        Exg_rate_disposed_system,Exg_q_rate_out_hospital,Exg_f_rate_out_hospital, ...
        Exg_q_rate_loss_furnace, Exg_d_rate_fan,Exg_d_rate_heatingwhole,Exg_d_rate_hospital);
    
    %% Table - State Properties
    [State_Table, states, Specific_Energy, Carried_Energy_Rate, Specific_Flow_Exergy, Flow_Exergy_Rate] = States_Fn(m_dot, cp,R, T0, T1, T2, T3, T4, P0, P1, P2, P3, P4);
    
    %% Table - Sub Systems (Components) Energy and Performance Rate
    [SubSystems_Energy_Table, SubSystems_Performance_Table]=SubSystems_EnergyRate_Fn(Energy_w_rate_in_fan,Energy_w_rate_in_heatingwhole,...
    Energy_w_rate_in_hospital, Energy_w_rate_out_fan, Energy_w_rate_out_heatingwhole,Energy_w_rate_out_hospital,...
    Energy_q_rate_in_fan, Energy_q_rate_in_heatingwhole, Energy_q_rate_in_hospital,...
    Energy_q_rate_out_fan, Energy_q_rate_out_heatingwhole, Energy_q_rate_out_hospital,...
    Energy_q_rate_loss_fan, Energy_q_rate_loss_heatingwhole, Energy_q_rate_loss_hospital,...
    ...
    Energy_rate_in_fan, Energy_rate_in_heatingwhole,Energy_rate_in_hospital,...
    Energy_utilized_rate_fan, Energy_utilized_rate_heatingwhole,Energy_utilized_rate_hospital,...
    ...
    Exg_rate_in_fan, Exg_rate_in_heatingwhole,Exg_rate_in_hospital,...
    Exg_utilized_rate_fan, Exg_utilized_rate_heatingwhole,Exg_utilized_rate_hospital,...
    Exg_q_rate_loss_fan, Exg_q_rate_loss_heatingwhole, Exg_q_rate_loss_hospital,...
    ...
    sigma_dot_p_fan, sigma_dot_p_heatingwhole, sigma_dot_p_hospital, Exg_d_rate_fan, Exg_d_rate_heatingwhole, Exg_d_rate_hospital,...
    eff_energy_fan, eff_energy_heatingwhole, eff_energy_hospital,eff_exergy_fan, eff_exergy_heatingwhole, eff_exergy_hospital);
    
    %% (For period) Total Energy/ Exergy/ Entropy production ** Needs Modification
    [SubSystems_Energy_Period_Table, SubSystems_Performance_Period_Table,Energy_in_fan_period,Energy_in_heatingwhole_period]=...
    SubSystems_EnergyPeriod_Fn(dt_analysis, Energy_w_rate_in_fan,Energy_w_rate_in_heatingwhole,...
    Energy_w_rate_in_hospital, Energy_w_rate_out_fan, Energy_w_rate_out_heatingwhole,Energy_w_rate_out_hospital,...
    Energy_q_rate_in_fan, Energy_q_rate_in_heatingwhole, Energy_q_rate_in_hospital,...
    Energy_q_rate_out_fan, Energy_q_rate_out_heatingwhole, Energy_q_rate_out_hospital,...
    Energy_q_rate_loss_fan, Energy_q_rate_loss_heatingwhole, Energy_q_rate_loss_hospital,...
    ...
    Energy_rate_in_fan, Energy_rate_in_heatingwhole,Energy_rate_in_hospital,...
    Energy_utilized_rate_fan, Energy_utilized_rate_heatingwhole,Energy_utilized_rate_hospital,...
    ...
    Exg_w_rate_in_fan,Exg_w_rate_in_heatingwhole,...
    Exg_w_rate_in_hospital, Exg_w_rate_out_fan, Exg_w_rate_out_heatingwhole,Exg_w_rate_out_hospital,...
    Exg_q_rate_in_fan, Exg_q_rate_in_heatingwhole, Exg_q_rate_in_hospital,...
    Exg_q_rate_out_fan, Exg_q_rate_out_heatingwhole, Exg_q_rate_out_hospital, ...
    Exg_q_rate_loss_fan, Exg_q_rate_loss_heatingwhole, Exg_q_rate_loss_hospital,...
    ...
    Exg_rate_in_fan, Exg_rate_in_heatingwhole,Exg_rate_in_hospital,...
    Exg_utilized_rate_fan, Exg_utilized_rate_heatingwhole,Exg_utilized_rate_hospital,...
    ...
    sigma_dot_p_fan, sigma_dot_p_heatingwhole, sigma_dot_p_hospital, Exg_d_rate_fan, Exg_d_rate_heatingwhole, Exg_d_rate_hospital,...
    eff_energy_fan, eff_energy_heatingwhole, eff_energy_hospital,eff_exergy_fan, eff_exergy_heatingwhole, eff_exergy_hospital);
    %% Energy Cost Per Period 
    [specific_volume_ng, natural_gas_energy_MJ_per_m3, Energy_in_fan_MJ_period, Energy_in_fan_kWh_period, electricity_cost_fan_period,...
        Energy_in_heatingwhole_MJ_period, natural_gas_volume_m3_period, fuel_cost_furnace_period]=Energy_Cost_Period_Fn(density_ng,...
        HHV_ng, Energy_in_fan_period,Energy_in_heatingwhole_period,electricity_cost_per_kWh,natural_gas_cost_per_m3);
    %% Sotoed in list Until Monthly Period is reached
    Energy_in_fan_MJ_period_list(jj) = Energy_in_fan_MJ_period;
    Energy_in_fan_kWh_period_list(jj) = Energy_in_fan_kWh_period;
    electricity_cost_fan_period_list(jj) = electricity_cost_fan_period;
    
    Energy_in_heatingwhole_MJ_period_list(jj) = Energy_in_heatingwhole_MJ_period;
    natural_gas_volume_m3_period_list(jj) = natural_gas_volume_m3_period;
    fuel_cost_furnace_period_list(jj)=fuel_cost_furnace_period;
        end
    %% End of looping over a month
    %% At the end of looping over a month, Monthly data is calculated to proceed monthly total cost
    %% Including Operator/ Capital and Energy Cost
     
        if ~strcmpi(analysis_period,'Winter')
            %% Monthly Energy Total
            Energy_in_fan_MJ_monthly(ii) = sum(Energy_in_fan_MJ_period_list);
            Energy_in_fan_kWh_monthly(ii) = sum(Energy_in_fan_kWh_period_list);
            electricity_cost_fan_monthly(ii) = sum(electricity_cost_fan_period_list);
            
            Energy_in_heatingwhole_MJ_monthly(ii) =sum(Energy_in_heatingwhole_MJ_period_list);
            natural_gas_volume_m3_monthly(ii) = sum(natural_gas_volume_m3_period_list);
            fuel_cost_furnace_monthly(ii) = sum(fuel_cost_furnace_period_list);
    
            %% Monthly Energy Cost Total
            months_period_analyzed = 1; %Monthly
            name_period_analyzed = [months_design{ii}  ' (Monthly)'];
            electricity_cost_fan_period_analyzed = electricity_cost_fan_monthly(ii);
            fuel_cost_furnace_period_analyzed = fuel_cost_furnace_monthly(ii);
    
            [operator_cost_period_monthly, total_cost_woCapital_period_monthly, additional_payment_CHEHP_period_monthly,additional_payment_regHE_period_monthly, ...
    total_cost_withCapital_period_monthly]=All_Cost_MonthlyOrWinter_CHEHP_Fn(name_period_analyzed,months_period_analyzed,operator_full_cost_per_year,operator_allocation,...
    electricity_cost_fan_period_analyzed, fuel_cost_furnace_period_analyzed,...
    capital_cost_fan, capital_cost_furnace, capital_cost_he, captial_cost_hp, capital_cost_regHE, annual_interest_rate,financing_period_years);
    
            %% Monthly Carbon Emissions Total
            Energy_in_fan_kWh_period_analyzed = Energy_in_fan_kWh_monthly(ii);
            natural_gas_volume_m3_period_analyzed = natural_gas_volume_m3_monthly(ii);
    
            [electricity_geneartaion_kWh_period_analyzed,electricity_fossil_CO2_emissions_period_analyzed,...
            total_CO2_emissions_period_analyzed,Carbon_Emissions_Table_Period]=CarbonEmissions_MonthlyOrWinter_Fn(name_period_analyzed,Mol_CO2,Mol_CH4,HHV_ng,eff_power_plant_ng,...
            density_ng,renewable_share_electricity,transmission_loss_electricity, ...
            Energy_in_fan_kWh_period_analyzed, natural_gas_volume_m3_period_analyzed);
            end
    
    end
    
    %% End of Over the season or Design months
    % Total Energy for the Winter - Monthly to Winter  // Or Winter Anlaysis
    if strcmpi(analysis_period,'Winter')
            Energy_in_fan_MJ_winter = sum(Energy_in_fan_MJ_period_list); %of course no need to sum - only one in list
            Energy_in_fan_kWh_winter = sum(Energy_in_fan_kWh_period_list);
            electricity_cost_fan_winter = sum(electricity_cost_fan_period_list);
            
            Energy_in_heatingwhole_MJ_winter =sum(Energy_in_heatingwhole_MJ_period_list);
            natural_gas_volume_m3_winter = sum(natural_gas_volume_m3_period_list);
            fuel_cost_furnace_winter = sum(fuel_cost_furnace_period_list);
    else
        Energy_in_fan_MJ_winter = sum(Energy_in_fan_MJ_monthly);
        Energy_in_fan_kWh_winter = sum(Energy_in_fan_kWh_monthly);
        electricity_cost_fan_winter = sum(electricity_cost_fan_monthly);
                
        Energy_in_heatingwhole_MJ_winter =sum(Energy_in_heatingwhole_MJ_monthly);
        natural_gas_volume_m3_winter = sum(natural_gas_volume_m3_monthly);
        fuel_cost_furnace_winter = sum(fuel_cost_furnace_monthly);
    end
    %% Total Energy Cost for the whole season - Winter
    months_period_analyzed = length(months_design);
    name_period_analyzed = 'Winter';
    electricity_cost_fan_period_analyzed = electricity_cost_fan_winter;
    fuel_cost_furnace_period_analyzed = fuel_cost_furnace_winter;
    
    [Cost_Table,operator_cost_period_winter, total_cost_woCapital_period_winter, additional_payment_CHEHP_period_winter, additional_payment_regHE_period_winter,...
    total_cost_withCapital_period_winter]=All_Cost_MonthlyOrWinter_CHEHP_Fn(name_period_analyzed,months_period_analyzed,operator_full_cost_per_year,operator_allocation,...
    electricity_cost_fan_period_analyzed, fuel_cost_furnace_period_analyzed,...
    capital_cost_fan, capital_cost_furnace, capital_cost_he, captial_cost_hp, capital_cost_regHE, annual_interest_rate,financing_period_years);

   
    %% Total Carbon Emissions for the whole season - Winter
    Energy_in_fan_kWh_period_analyzed = Energy_in_fan_kWh_winter;
    natural_gas_volume_m3_period_analyzed = natural_gas_volume_m3_winter;
            
    [electricity_geneartaion_kWh_winter_winter,electricity_fossil_CO2_emissions_winter,...
        natural_gas_CO2_emissions_winter,total_CO2_emissions_winter,Carbon_Emissions_Table_winter]=...
        CarbonEmissions_MonthlyOrWinter_Fn(name_period_analyzed,Mol_CO2,Mol_CH4,HHV_ng,eff_power_plant_ng,...
    density_ng,renewable_share_electricity,transmission_loss_electricity, ...
    Energy_in_fan_kWh_period_analyzed, natural_gas_volume_m3_period_analyzed); 

    %% For temperature change effects >> Compact variables for plot

    %%  *** Assumptions - Energy and Exergy Efficiency do not change throughout the system
    
    eff_energy_system_winter_list_Teffect(tt) = eff_energy_system*100; %Add
    eff_exergy_system_winter_list_Teffect(tt) = eff_exergy_system*100; %Add

    % Fan Power and Fuel Cost
    Energy_in_fan_MJ_winter_list_Teffect(tt) = Energy_in_fan_MJ_winter; %Add
    Energy_in_fan_kWh_winter_list_Teffect(tt) = Energy_in_fan_kWh_winter;
    electricity_cost_fan_winter_list_Teffect(tt) = electricity_cost_fan_winter;

    Energy_in_heatingwhole_MJ_winter_list_Teffect(tt)=Energy_in_heatingwhole_MJ_winter; %Add
    natural_gas_volume_m3_winter_list_Teffect(tt) = natural_gas_volume_m3_winter;
    fuel_cost_furnace_winter_list_Teffect(tt) = fuel_cost_furnace_winter;

    total_cost_withCapital_period_winter_list_Teffect(tt) = total_cost_withCapital_period_winter;
    
    % CO2 Emissions
    electricity_fossil_CO2_emissions_winter_list_Teffect(tt) = electricity_fossil_CO2_emissions_winter;
    natural_gas_CO2_emissions_winter_list_Teffect(tt)=natural_gas_CO2_emissions_winter;
    
    total_CO2_emissions_winter_list_Teffect(tt) = total_CO2_emissions_winter;

    % Exergy destruction and Exe
    Energetic_efficiency_winter_list_Teffect(tt,:) = SubSystems_Performance_Period_Table{:,10}'; % 7th Column is energetic efficiency
    Exergetic_efficiency_winter_list_Teffect(tt,:) = SubSystems_Performance_Period_Table{:,11}'; % 7th Column is energetic efficiency

end
%% Add Next for Overall Efficiency
    Energetic_efficiency_winter_list_Teffect = [Energetic_efficiency_winter_list_Teffect eff_energy_system_winter_list_Teffect']; %Add
    Exergetic_efficiency_winter_list_Teffect= [Exergetic_efficiency_winter_list_Teffect eff_exergy_system_winter_list_Teffect']; %Add


%% Saving Variables
T0_list_option1 = T0_list;
save('HeatingSystemData_option1.mat', ...
    'T0_list_option1', ...
    'Energy_in_heatingwhole_MJ_winter_list_Teffect', ...
    'fuel_cost_furnace_winter_list_Teffect', ...
    'total_cost_withCapital_period_winter_list_Teffect', ...
    'natural_gas_CO2_emissions_winter_list_Teffect', ...
    'Energetic_efficiency_winter_list_Teffect', ...
    'Exergetic_efficiency_winter_list_Teffect');

%% Excel Results
Eff_Labels = {'The Whole System'};
Eff_Table = table(Eff_Labels, eff_energy_system*100, eff_exergy_system*100, 'VariableNames', {'System', 'Energy Efficiency/EUF [%]','Exergy Efficiency [%]'});

% Define the file path for the Excel output
filePath = 'Hospital_Project_Tables_Option_1.xlsx'; % Change this to your desired path
% Define the sheet names for each table
sheetNames = {'State_Table', ...
               'exergy_sources_table', ...
              'exergy_disposition_table', ...
              'SubSystems_Energy_Table', ...
              'SubSystems_Performance_Table', ...
              'SubSystems_Energy_Period_Table', ...
              'SubSystems_Perform_Period_Table', ...
              'The Whole System',...
              'Cost_Table', ...
              'Carbon_Emissions_Table_winter'};

% Define the table variables
tables = {State_Table, ...
    exergy_sources_table, ...
          exergy_disposition_table, ...
          SubSystems_Energy_Table, ...
          SubSystems_Performance_Table, ...
          SubSystems_Energy_Period_Table, ...
          SubSystems_Performance_Period_Table, ...
          Eff_Table,...
          Cost_Table, ...
          Carbon_Emissions_Table_winter};

% Write each table to the specified sheet
for i = 1:length(tables)
    writetable(tables{i}, filePath, 'Sheet', sheetNames{i});
end

disp('All tables have been written to the specified Excel file.');

%% For GUI

% Example plot data for grouped bar chart
    plotData.x = T0_list; % X-axis (ambient temperature)
    plotData.y = Exergetic_efficiency_winter_list_Teffect(:,1:4); % Replace with actual efficiency data
    plotData.legend = {'Fan','Furnace (the whole heating process by Furnace)', 'Hospital Distribution','Overall System (Fan + ALL)'};

end