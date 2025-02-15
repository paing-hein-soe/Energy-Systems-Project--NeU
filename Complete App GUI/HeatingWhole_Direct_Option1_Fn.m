function [Energy_w_rate_in_heatingwhole, Energy_q_rate_in_heatingwhole, Energy_rate_in_heatingwhole,...
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
    Exg_utilized_rate_heatingair, Exg_q_rate_loss_heatingair, Exg_q_rate_loss_furnace)
%% Whole System of Furnace Heating - Heating Air + Furnace

% Energy Analysis 
Energy_q_rate_in_heatingwhole = Energy_q_rate_in_furnace;
Energy_w_rate_in_heatingwhole=Energy_w_rate_in_furnace;
Energy_rate_in_heatingwhole = Energy_q_rate_in_heatingwhole+Energy_w_rate_in_heatingwhole; 

Energy_q_rate_out_heatingwhole= Energy_q_rate_out_heatingair;
Energy_w_rate_out_heatingwhole = Energy_w_rate_out_heatingair;
Energy_rate_out_heatingwhole = Energy_q_rate_out_heatingwhole+Energy_w_rate_out_heatingwhole;

Energy_utilized_rate_heatingwhole=Energy_utilized_rate_heatingair; %Only heating air is utilized
Energy_q_rate_loss_heatingwhole = Energy_q_rate_loss_heatingair + Energy_q_rate_loss_furnace; 

%% Entropy production
sigma_dot_p_heatingwhole = sigma_dot_p_furnace+sigma_dot_p_heatingair;
Exg_d_rate_heatingwhole = Exg_d_rate_furnace + Exg_d_rate_heatingair;

%% Exergy Flow from 2 -3 Whole Heating Process
exg_f_in_heatingwhole=exg_f_in_heatingair;
exg_f_out_heatingwhole=exg_f_out_heatingair;
Exg_f_rate_in_heatingwhole = Exg_f_rate_in_heatingair;
Exg_f_rate_out_heatingwhole = Exg_f_rate_out_heatingair;

% Exergy in/out
Exg_q_rate_in_heatingwhole = Exg_q_rate_in_furnace;
Exg_w_rate_in_heatingwhole =Exg_w_rate_in_furnace;
Exg_rate_in_heatingwhole=Exg_q_rate_in_heatingwhole+Exg_w_rate_in_heatingwhole;

Exg_q_rate_out_heatingwhole = Exg_q_rate_out_heatingair;
Exg_w_rate_out_heatingwhole = Exg_w_rate_out_heatingair;
Exg_rate_out_heatingwhole=Exg_q_rate_out_heatingwhole+Exg_w_rate_out_heatingwhole;
Exg_q_rate_loss_heatingwhole = Exg_q_rate_loss_heatingair + Exg_q_rate_loss_furnace;

Exg_d_rate_heatingwhole_exgbalance= (Exg_q_rate_in_heatingwhole-Exg_q_rate_out_heatingwhole) +(Exg_w_rate_in_heatingwhole-Exg_w_rate_out_heatingwhole)...
    +(-Exg_q_rate_loss_heatingwhole)+(Exg_f_rate_in_heatingwhole -Exg_f_rate_out_heatingwhole);

Exg_utilized_rate_heatingwhole=Exg_utilized_rate_heatingair;

%% Efficiency
eff_energy_heatingwhole= Energy_utilized_rate_heatingwhole/Energy_rate_in_heatingwhole;
eff_exergy_heatingwhole= Exg_utilized_rate_heatingwhole/Exg_rate_in_heatingwhole;

%% Energy/Exergy Results -Heating
fprintf('\nThe Whole Heating System - Furnace and Heating Air- Results:\n');

fprintf('Exergy Destroyed Rate in the whole Heating System (W): %.2f\n', Exg_d_rate_heatingwhole);
fprintf('Exergy Destroyed Rate (Check- Exergy Balance) in the whole Heating System (W): %.2f\n', Exg_d_rate_heatingwhole_exgbalance);

end