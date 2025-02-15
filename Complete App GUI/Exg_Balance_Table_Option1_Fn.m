function [Exg_rate_in_system_check,Exg_rate_disposed_system_check,exergy_sources_table, exergy_disposition_table]=...
    Exg_Balance_Table_Option1_Fn(m_dot,Exg_rate_in_system,Exg_f_rate_in_fan,Exg_w_rate_in_fan,Exg_q_rate_in_heatingwhole,...
    Exg_rate_disposed_system,Exg_q_rate_out_hospital,Exg_f_rate_out_hospital, ...
    Exg_q_rate_loss_furnace,Exg_d_rate_fan,Exg_d_rate_heatingwhole,Exg_d_rate_hospital)
  

    %% Extension for Exergy Balance Table
    % Data for Sources of Exergy
    Exg_rate_in_system_check = Exg_f_rate_in_fan+Exg_w_rate_in_fan+Exg_q_rate_in_heatingwhole;
    sources_exergy = {
  
        'Exegy of Flow to Fan', Exg_f_rate_in_fan, 'W', Exg_f_rate_in_fan/m_dot, 'kJ/kg';
        'Exergy of Work Input to Fan', Exg_w_rate_in_fan, 'W', Exg_w_rate_in_fan/m_dot, 'kJ/kg';
        'Exergy of Heat Input to Furnace', Exg_q_rate_in_heatingwhole, 'W',Exg_q_rate_in_heatingwhole/m_dot, 'kJ/kg';
        'Total Exergy Rate Input', Exg_rate_in_system_check,'W',Exg_rate_in_system_check/m_dot,'kJ/kg';
        'Total Exergy Rate Input(Check)',Exg_rate_in_system,'W',Exg_rate_in_system/m_dot,'kJ/kg'
    };
    
    % Data for Exergy Disposition
    Exg_rate_disposed_system_check = Exg_q_rate_out_hospital+Exg_q_rate_loss_furnace+Exg_f_rate_out_hospital+...
        Exg_d_rate_fan+ Exg_d_rate_heatingwhole+Exg_d_rate_hospital;
    disposition_exergy = {
        
        'Useful Exergy by Hospital Distribution', Exg_q_rate_out_hospital, 'W', Exg_q_rate_out_hospital/m_dot, 'kJ/kg';
        'Exergy of Heat Loss by Furnace', Exg_q_rate_loss_furnace, 'W', Exg_q_rate_loss_furnace/m_dot, 'kJ/kg';
        'Exergy of Flow Out at Hospital Exit', Exg_f_rate_out_hospital, 'W', Exg_f_rate_out_hospital/m_dot, 'kJ/kg';
        'Exergy Destruction (Fan)', Exg_d_rate_fan, 'W', Exg_d_rate_fan/m_dot, 'kJ/kg';
        'Exergy Destruction (Furnace-Heating Whole)', Exg_d_rate_heatingwhole, 'W', Exg_d_rate_heatingwhole/m_dot, 'kJ/kg';
        'Exergy Destruction (Hospital)', Exg_d_rate_hospital, 'W', Exg_d_rate_hospital/m_dot, 'kJ/kg';
        'Total Exergy Rate Disposed', Exg_rate_disposed_system_check, 'W', Exg_rate_disposed_system_check/m_dot, 'kJ/kg';
        'Total Exergy Rate Disposed(Check)', Exg_rate_disposed_system, 'W', Exg_rate_disposed_system/m_dot, 'kJ/kg';

    };
    
    % Convert data to MATLAB table
    exergy_sources_table = cell2table(sources_exergy, ...
        'VariableNames', {'Sources of Exergy','Exergy Rate', 'Unit1','Specific Exergy','Unit2'});
    
    exergy_disposition_table = cell2table(disposition_exergy, ...
        'VariableNames', {'Exergy Disposition','Exergy Rate', 'Unit1','Specific Exergy','Unit2'});
    
    % Display tables
    fprintf('Table A1: Exergy Balance (Sources of Exergy)\n');
    disp(exergy_sources_table);
    
    fprintf('Table A2: Exergy Balance (Exergy Disposition)\n');
    disp(exergy_disposition_table);

end