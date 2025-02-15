function [Energy_rate_in_system, Exg_rate_disposed_system, Energy_utilized_rate_system,...
    Exg_rate_in_system, Exg_utilized_rate_system, Exg_d_rate_system, Exg_q_rate_loss_system,...
    eff_energy_system, eff_exergy_system]= Overall_EneExeEff_Option4_Fn(...
    Energy_rate_in_fan, Energy_rate_in_heatingwhole,Energy_utilized_rate_hospital,...
    Exg_rate_in_fan, Exg_rate_in_heatingwhole,Exg_utilized_rate_hospital,...
    Exg_d_rate_fan, Exg_d_rate_heatingwhole, Exg_d_rate_hospital,...
    Exg_q_rate_loss_fan, Exg_q_rate_loss_heatingwhole, Exg_q_rate_loss_hospital, Exg_f_rate_out_regHE)

%% Total Energy Rate Utilized
Energy_utilized_rate_system = Energy_utilized_rate_hospital;
Energy_rate_in_system = Energy_rate_in_fan + Energy_rate_in_heatingwhole;

Exg_utilized_rate_system = Exg_utilized_rate_hospital; %Exg_f_rate_out_furnace- Exg_f_rate_out_hospital; % Utilized exergy (W)
Exg_rate_in_system = Exg_rate_in_fan + Exg_rate_in_heatingwhole;
Exg_d_rate_system = Exg_d_rate_fan + Exg_d_rate_heatingwhole + Exg_d_rate_hospital; 
Exg_q_rate_loss_system = Exg_q_rate_loss_fan + Exg_q_rate_loss_heatingwhole + Exg_q_rate_loss_hospital; 

% Disposed exergy
Exg_rate_disposed_system = Exg_q_rate_loss_system+Exg_d_rate_system+(Exg_utilized_rate_system-Exg_d_rate_hospital)+ Exg_f_rate_out_regHE; %% RegHE
%% Energy/Exergy Efficiency
eff_energy_system = Energy_utilized_rate_system / Energy_rate_in_system;
eff_exergy_system = Exg_utilized_rate_system / Exg_rate_in_system;

%% Results- Overall System Efficiencies
fprintf('\nOverall System Efficiency Results:');
fprintf('Overall First-Law (Energy) Efficiency: %.2f%%\n', eff_energy_system * 100);
fprintf('Overall Second-Law (Exergy) Efficiency: %.2f%%\n', eff_exergy_system * 100);
end