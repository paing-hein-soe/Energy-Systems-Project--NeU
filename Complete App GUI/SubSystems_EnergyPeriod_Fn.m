function [SubSystems_Energy_Period_Table, SubSystems_Performance_Period_Table,Energy_in_fan_period,Energy_in_heatingwhole_period]=...
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
    eff_energy_fan, eff_energy_heatingwhole, eff_energy_hospital,eff_exergy_fan, eff_exergy_heatingwhole, eff_exergy_hospital)

    %% 
    % Fan Analysis (For the period input - eg. Winter or November)
    
    Energy_w_in_fan_period = Energy_w_rate_in_fan * dt_analysis; % Exergy input (J/year)
    Energy_w_in_heatingwhole_period = Energy_w_rate_in_heatingwhole * dt_analysis; % Exergy input (J/year)
    Energy_w_in_hospital_period= Energy_w_rate_in_hospital * dt_analysis; % Exergy input (J/year)
    
    Energy_w_out_fan_period = Energy_w_rate_out_fan * dt_analysis; % Exergy input (J/year)
    Energy_w_out_heatingwhole_period = Energy_w_rate_out_heatingwhole * dt_analysis; % Exergy input (J/year)
    Energy_w_out_hospital_period= Energy_w_rate_out_hospital * dt_analysis; % Exergy input (J/year)
    
    Energy_q_in_fan_period = Energy_q_rate_in_fan * dt_analysis; % Exergy input (J/year)
    Energy_q_in_heatingwhole_period = Energy_q_rate_in_heatingwhole * dt_analysis; % Exergy input (J/year)
    Energy_q_in_hospital_period= Energy_q_rate_in_hospital * dt_analysis; % Exergy input (J/year)
    
    Energy_q_out_fan_period = Energy_q_rate_out_fan * dt_analysis; % Exergy input (J/year)
    Energy_q_out_heatingwhole_period = Energy_q_rate_out_heatingwhole * dt_analysis; % Exergy input (J/year)
    Energy_q_out_hospital_period= Energy_q_rate_out_hospital * dt_analysis; % Exergy input (J/year)

    Energy_q_loss_fan_period = Energy_q_rate_loss_fan * dt_analysis; % Exergy input (J/year)
    Energy_q_loss_heatingwhole_period = Energy_q_rate_loss_heatingwhole * dt_analysis; % Exergy input (J/year)
    Energy_q_loss_hospital_period= Energy_q_rate_loss_hospital * dt_analysis; % Exergy input (J/year)

%%
     %% Table: Energy Flows Related to Components
     SubSystems_Label = {'Fan'; 'The Whole Heating System by Furnace'; 'Hospital Air Distribution'};
     Power_In_SubSystems = [Energy_w_in_fan_period; Energy_w_in_heatingwhole_period; Energy_w_in_hospital_period];
     Power_Out_SubSystems = [Energy_w_out_fan_period; Energy_w_out_heatingwhole_period; Energy_w_out_hospital_period];
     Heat_In_SubSystems = [Energy_q_in_fan_period; Energy_q_in_heatingwhole_period; Energy_q_in_hospital_period];
     Heat_Out_SubSystems = [Energy_q_out_fan_period; Energy_q_out_heatingwhole_period; Energy_q_out_hospital_period]; % - Q_dot_ for heating will change
     Heat_Loss_SubSystems = [Energy_q_loss_fan_period; Energy_q_loss_heatingwhole_period; Energy_q_loss_hospital_period]; % - Q_dot_ for heating will change
     SubSystems_Energy_Period_Table = table(SubSystems_Label, Power_In_SubSystems, Power_Out_SubSystems, ...
         Heat_In_SubSystems, Heat_Out_SubSystems,Heat_Loss_SubSystems,...
            'VariableNames', {'Subsystem (Components)','Power In [J/period]', 'Power Out [J/period]',...
            'Heat In [J/period]','Heat Out [J/period]','Heat Loss [J/period]'});
        
     disp('Table 4: Energy (Total for Specified Period) Related to Components (or Sub-systems)');
     disp(SubSystems_Energy_Period_Table);      
        

%% Exergy  Each component

    sigma_p_fan_period = sigma_dot_p_fan * dt_analysis; % Entropy production (J/K per period)
    sigma_p_heatingwhole_period = sigma_dot_p_heatingwhole * dt_analysis; % Entropy production (J/K per period)
    sigma_p_hospital_period = sigma_dot_p_hospital * dt_analysis; % Entropy production (J/K per period)

    Exg_d_fan_period = Exg_d_rate_fan * dt_analysis; % Exergy Destruction (J per period)
    Exg_d_heatingwhole_period = Exg_d_rate_heatingwhole * dt_analysis; % Exergy Destruction (J per period)
    Exg_d_hospital_period = Exg_d_rate_hospital * dt_analysis; % Exergy Destruction (J per period)

    Energy_in_fan_period = Energy_rate_in_fan * dt_analysis; % Energy IN (J per period)
    Energy_in_heatingwhole_period = Energy_rate_in_heatingwhole * dt_analysis; % Energy IN (J per period)
    Energy_in_hospital_period = Energy_rate_in_hospital * dt_analysis; % Energy IN (J per period)

    Energy_utilized_fan_period = Energy_utilized_rate_fan * dt_analysis; % Energy Utilized (J per period)
    Energy_utilized_heatingwhole_period = Energy_utilized_rate_heatingwhole * dt_analysis; % Energy (J per period)
    Energy_utilized_hospital_period = Energy_utilized_rate_hospital * dt_analysis; % Energy Utilized (J per period)
    
    Exg_in_fan_period = Exg_rate_in_fan * dt_analysis; % Energy IN (J per period)
    Exg_in_heatingwhole_period = Exg_rate_in_heatingwhole * dt_analysis; % Energy IN (J per period)
    Exg_in_hospital_period = Exg_rate_in_hospital * dt_analysis; % Energy IN (J per period)

    Exg_loss_fan_period = Exg_q_rate_loss_fan * dt_analysis; % Energy IN (J per period)
    Exg_loss_heatingwhole_period = Exg_q_rate_loss_heatingwhole * dt_analysis; % Energy IN (J per period)
    Exg_loss_hospital_period = Exg_q_rate_loss_hospital * dt_analysis; % Energy IN (J per period)


    Exg_utilized_fan_period = Exg_utilized_rate_fan * dt_analysis; % Energy Utilized (J per period)
    Exg_utilized_heatingwhole_period = Exg_utilized_rate_heatingwhole * dt_analysis; % Energy (J per period)
    Exg_utilized_hospital_period = Exg_utilized_rate_hospital * dt_analysis; % Energy Utilized (J per period)

    %% For Table 
    sigma_p_period = [sigma_p_fan_period; sigma_p_heatingwhole_period; sigma_p_hospital_period];
    Exg_d_period = [Exg_d_fan_period; Exg_d_heatingwhole_period; Exg_d_hospital_period];

    Energy_in_period = [Energy_in_fan_period; Energy_in_heatingwhole_period; Energy_in_hospital_period];
    Energy_utilized_period = [Energy_utilized_fan_period; Energy_utilized_heatingwhole_period;Energy_utilized_hospital_period];
    Energy_loss_period=Heat_Loss_SubSystems;

    Exg_in_period = [Exg_in_fan_period; Exg_in_heatingwhole_period;Exg_in_hospital_period];
    Exg_utilized_period = [Exg_utilized_fan_period; Exg_utilized_heatingwhole_period;Exg_utilized_hospital_period];
    Exg_loss_period = [Exg_loss_fan_period; Exg_loss_heatingwhole_period;Exg_loss_hospital_period];

    Energy_utilization_percent_period = [eff_energy_fan * 100; eff_energy_heatingwhole * 100; eff_energy_hospital*100];
    Exg_utilization_percent_period = [eff_exergy_fan * 100; eff_exergy_heatingwhole * 100; eff_exergy_hospital*100];
        
    %% Table: Performance Related to Components

    SubSystems_Performance_Period_Table = table(SubSystems_Label, sigma_p_period, Exg_d_period, ...
            Energy_in_period, Energy_utilized_period, Energy_loss_period, Exg_in_period, Exg_utilized_period, ...
            Energy_loss_period, Energy_utilization_percent_period, Exg_utilization_percent_period,...
            'VariableNames', {'Subsystem (Components)','Entropy Production [J/K-period]','Exergy Destruction [J/period]',...
            'Energy Input [J/period]','Energy Utilized [J/period]', 'Energy Loss [J/period]','Exergy Input [J/period]','Exergy Utilized [J/period]'...
            'Exergy Loss [J/period]','Energy Utilization Efficiency [%]','Exergetic Efficiency [%]'});
        
        disp('Table 5: Performance (Total for Specified Period) Related to Components (or Sub-systems)');
        disp(SubSystems_Performance_Period_Table);        

end
