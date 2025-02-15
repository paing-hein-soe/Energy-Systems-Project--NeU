function [Energy_w_rate_in_furnace, Energy_q_rate_in_furnace, Energy_rate_in_furnace, ...
    Energy_w_rate_out_furnace, Energy_q_rate_out_furnace, Energy_rate_out_furnace, Energy_q_rate_loss_furnace,...
    ...
    Exg_w_rate_in_furnace, Exg_q_rate_in_furnace, Exg_rate_in_furnace, ...
    Exg_w_rate_out_furnace, Exg_q_rate_out_furnace, Exg_rate_out_furnace, Exg_q_rate_loss_furnace, ...
    ...
    ...
    sigma_dot_p_furnace, Exg_d_rate_furnace, Exg_d_rate_furnace_exgbalance, ...
    ...
    Energy_utilized_rate_furnace, Exg_utilized_rate_furnace, eff_exergy_furnace, eff_energy_furnace]=Furnace_Direct_Option1_Fn...
(T0, T_flame, eta_furnace, Energy_q_rate_in_heatingair, Energy_w_rate_in_heatingair, Exg_q_rate_in_heatingair, Exg_w_rate_in_heatingair)


%% Heat Input from Fuel for Furnace
Q_dot_fuel = Energy_q_rate_in_heatingair/ eta_furnace; % Heat input from fuel (W)
Q_dot_exhaustfuel = Q_dot_fuel - Energy_q_rate_in_heatingair;

%% Energy Analysis for Furnace 
Energy_q_rate_in_furnace = Q_dot_fuel;
Energy_w_rate_in_furnace = 0;
Energy_rate_in_furnace = Energy_q_rate_in_furnace + Energy_w_rate_in_furnace;
 
Energy_q_rate_out_furnace = Energy_q_rate_in_heatingair; %Linked to Heating Air Sytem
Energy_w_rate_out_furnace = Energy_w_rate_in_heatingair;
Energy_rate_out_furnace = Energy_q_rate_out_furnace + Energy_w_rate_out_furnace;

Energy_utilized_rate_furnace = Energy_rate_out_furnace;
Energy_q_rate_loss_furnace = Q_dot_exhaustfuel;

%% Entropy Production
% No Mass Flow ==> No Delta_S required
sigma_dot_p_furnace = -(Energy_q_rate_in_furnace/T_flame - Energy_q_rate_out_furnace/T_flame - Energy_q_rate_loss_furnace/T_flame);
Exg_d_rate_furnace = T0*sigma_dot_p_furnace;

%% Exergy Analysis for Furnace 

% No Exergy Flow : exg_f & Exg_f

% Exergy in/out
Exg_q_rate_in_furnace = (1-T0/T_flame)*Q_dot_fuel;
Exg_w_rate_in_furnace = 0 ;
Exg_rate_in_furnace = Exg_q_rate_in_furnace + Exg_w_rate_in_furnace;

Exg_q_rate_out_furnace = Exg_q_rate_in_heatingair;
Exg_w_rate_out_furnace = Exg_w_rate_in_heatingair;
Exg_rate_out_furnace = Exg_q_rate_out_furnace + Exg_w_rate_out_furnace; 
Exg_q_rate_loss_furnace = (1-T0/T_flame)* Q_dot_exhaustfuel;

Exg_d_rate_furnace_exgbalance= (Exg_q_rate_in_furnace-Exg_q_rate_out_furnace) +...
    +(-Exg_q_rate_loss_furnace)+(Exg_w_rate_in_furnace-Exg_w_rate_out_furnace);
 
Exg_utilized_rate_furnace = Exg_rate_out_furnace;
% Efficiency
eff_energy_furnace= Energy_utilized_rate_furnace/Energy_rate_in_furnace;
eff_exergy_furnace= Exg_utilized_rate_furnace/Exg_rate_in_furnace;

%% Display Results
fprintf('\nFurnace System - Direct Heating - Option 1- Results:\n');
fprintf('Heat rate input to Furance by burning fuel(W): %.2f\n', Q_dot_fuel);

fprintf('Exergy Destroyed Rate in Furnace (W): %.2f\n', Exg_d_rate_furnace);
fprintf('Exergy Destroyed Rate in Furnace (W) (checked by Exergy Balnace): %.2f\n', Exg_d_rate_furnace_exgbalance);

%fprintf('Option 1- Direct Heating Furnace System - First-Law (Energy) Efficiency: %.2f%%\n', eff_energy_furnace * 100);
%fprintf('Option 1- Direct Heating Furnace System - Second-Law (Exergy) Efficiency: %.2f%%\n', eff_exergy_furnace * 100);

end



