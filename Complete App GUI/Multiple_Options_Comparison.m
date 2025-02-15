clc; clearvars; close all;
%% Saving Variables
line_type = {'-', '-.', ':','--'};
line_color = {'r','m','b','g'};
marker_type = {'o','.','+','*'};

% Define variable names to load and plot dynamically
variable_names = {
    'T0_list_option1', ...
    'Energy_in_heatingwhole_MJ_winter_list_Teffect', ...
    'fuel_cost_furnace_winter_list_Teffect', ...
    'total_cost_withCapital_period_winter_list_Teffect', ...
    'natural_gas_CO2_emissions_winter_list_Teffect', ...
    'Energetic_efficiency_winter_list_Teffect(:,end)', ...
    'Exergetic_efficiency_winter_list_Teffect(:,end)'
};

% Define properties for plotting each variable
plot_titles = {
    'Ambient Temperature [°C]', ...
    'Fuel Energy Input [MJ]', ...
    'Fuel Cost', ...
    'Total Cost with Capital', ...
    'Natural Gas CO2 Emissions', ...
    'Overall Energy Efficiency [%]', ...
    'Overall Exergy Efficiency [%]'
};

y_labels = {
    'Temperature [°C]', ...
    'Energy [MJ]', ...
    'Cost [USD]', ...
    'Cost [USD]', ...
    'Emissions [kg]', ...
    'Efficiency [%]', ...
    'Efficiency [%]'
};
legend_title = {'Option 1: Direct Heating by Furnace only', ...
          'Option 2: Furnace and CHEHP', ...
          'Option 3: Furnace & CHEHP Heat Engine Recovery', ...
  sprintf(['Option 4: Furnace, CHEHP with heat recovery from Heat Engine \n' ...
                   'and Hospital Discharge'])};
% Loop over each variable to create separate figures
for jj = 2:length(variable_names) % Skip T0_list_option1 since it's an x-axis variable
    figure; % Create a new figure for each variable
    
    % Loop over the options (1 to 4)
    for ii = 1:4
        load(['HeatingSystemData_option' num2str(ii) '.mat']); % Load the data
        
        % Extract data dynamically using variable names
        x_data = eval(['T0_list_option' num2str(ii)]); % x-axis data
        y_data = eval(variable_names{jj});           % y-axis data
        
        % Plot the data with specified styles
        plot(x_data, y_data, ...
            [line_type{ii} line_color{ii} marker_type{ii}], ...
            'LineWidth', 1.5); hold on;
        if ii==1
            xticks(x_data); % Set x-axis ticks to match x_data;
        end
    end
    % Customize x-axis ticks
    
    % Add labels, title, and legend
    xlabel('Ambient Temperature [°C]');
    ylabel(y_labels{jj});
    title(['Comparison of ' plot_titles{jj} ' Across Options']);
    hlegend = legend(legend_title);
    set(hlegend, 'FontSize', 9);
    grid on;
end
%%
%% Additional Cost Analysis

ind = find(T0_list_option1==8); % 8 degree Celsius
kk=1;
load(['HeatingSystemData_option' num2str(kk) '.mat']); % Load the data
fuel_cost_option1 = fuel_cost_furnace_winter_list_Teffect(ind);

capital_cost_he = 100000; % $100,000 
capital_cost_hp = 50000; % $50,000
capital_cost_regHE = 10000; %$10,000
capital_cost_chehp = capital_cost_he + capital_cost_hp; % $100,000 
capital_cost_chehpRegHE = capital_cost_chehp+ capital_cost_regHE;


load('HeatingSystemData_option1.mat'); % Load the data
fuel_cost_option1 = fuel_cost_furnace_winter_list_Teffect(ind);
fuel_savings_option1 = 0;

load('HeatingSystemData_option2.mat'); % Load the data
fuel_cost_option2 = fuel_cost_furnace_winter_list_Teffect(ind);
fuel_savings_option2 =  fuel_cost_option1 -fuel_cost_option2;
savings_ratio_addcapital_option2 = fuel_savings_option2 / capital_cost_chehp;
breakeven_winterseasons_option2 = 1/savings_ratio_addcapital_option2;

load('HeatingSystemData_option3.mat'); % Load the data
fuel_cost_option3 = fuel_cost_furnace_winter_list_Teffect(ind);
fuel_savings_option3 = fuel_cost_option1-fuel_cost_option3;
savings_ratio_addcapital_option3 = fuel_savings_option3 / capital_cost_chehp;
breakeven_winterseasons_option3 = 1/savings_ratio_addcapital_option3;

load('HeatingSystemData_option4.mat'); % Load the data
fuel_cost_option4 = fuel_cost_furnace_winter_list_Teffect(ind);
fuel_savings_option4 = fuel_cost_option1-fuel_cost_option4;
savings_ratio_addcapital_option4 = fuel_savings_option4/ capital_cost_chehpRegHE;
breakeven_winterseasons_option4 = 1/savings_ratio_addcapital_option4;

breakeven_winterseasons = [breakeven_winterseasons_option2;...
    breakeven_winterseasons_option3;breakeven_winterseasons_option4];
figure; plot(breakeven_winterseasons)

figure; 


