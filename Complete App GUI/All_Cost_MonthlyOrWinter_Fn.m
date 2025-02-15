function [operator_cost_period_analyzed, total_cost_woCapital_period_analyzed, total_cost_withCapital_period_analyzed]=All_Cost_MonthlyOrWinter_Fn(name_period_analyzed,months_period_analyzed,operator_full_cost_per_year,operator_allocation,...
    electricity_cost_fan_period_analyzed, fuel_cost_furnace_period_analyzed,...
    capital_cost_fan, capital_cost_furnace, annual_interest_rate,financing_period_years)
%% Considering Operator Cost 
months_period = months_period_analyzed ; % Monthly analyzed = 1; Seasonal => 6

operator_full_cost_per_month = operator_full_cost_per_year/12; %$/month
operator_full_cost_period_analyzed = operator_full_cost_per_month * months_period;

%operator_full_cost_per_month = operator_full_cost_per_month * months_period; % $/month

% Step 3: Operator Cost
operator_cost_period_analyzed = operator_full_cost_period_analyzed * operator_allocation; % 20% of operator salary ($)

% Step 4: Total Seasonal Cost
total_cost_woCapital_period_analyzed = electricity_cost_fan_period_analyzed + fuel_cost_furnace_period_analyzed + operator_cost_period_analyzed; % Total operating cost ($)

% Display Results
%{
fprintf('All Cost Analysis Results in "%s"\n',name_period_analyzed);
fprintf('Electricity Cost for Fan: $%.2f\n', electricity_cost_fan_period_analyzed);
fprintf('Fuel Cost for Heating by Furnace: $%.2f\n', fuel_cost_furnace_period_analyzed);
fprintf('Operator Cost: $%.2f\n', operator_cost_period_analyzed);
fprintf('Total Cost (without Capital Cost) for the whole operation: $%.2f\n', total_cost_woCapital_period_analyzed);
%}

%% Rough Estimate for Capital Cost
% Annualized Cost Calculation (Using Simple Annual Payment Formula)

annual_payment_fan = capital_cost_fan * (annual_interest_rate * (1 + annual_interest_rate)^financing_period_years) / ...
                     ((1 + annual_interest_rate)^financing_period_years - 1);
annual_payment_furnace = capital_cost_furnace * (annual_interest_rate * (1 + annual_interest_rate)^financing_period_years) / ...
                         ((1 + annual_interest_rate)^financing_period_years - 1);

% Adjust for Monthly or Winter 
payment_fan_period_analyzed = annual_payment_fan * (months_period / 12); % Seasonal fan cost
payment_furnace_period_analyzed = annual_payment_furnace * (months_period / 12); % Seasonal furnace cost

% All Costs
total_cost_withCapital_period_analyzed = total_cost_woCapital_period_analyzed + payment_fan_period_analyzed + payment_furnace_period_analyzed;
%fprintf('Total Cost (with Capital Cost) for the whole operation: $%.2f\n', total_cost_withCapital_period_analyzed);

%% Table:
% Create a table for cost analysis
Cost_Labels = {
    'Electricity Cost for Fan';
    'Fuel Cost for Heating by Furnace';
    'Operator Cost';
    'Total Cost (without Capital Cost)';
    'Total Cost (with Capital Cost)'
};

Cost_Values = [
    electricity_cost_fan_period_analyzed;
    fuel_cost_furnace_period_analyzed;
    operator_cost_period_analyzed;
    total_cost_woCapital_period_analyzed;
    total_cost_withCapital_period_analyzed
];

% Combine into a table
Cost_Table = table(Cost_Labels, Cost_Values, ...
    'VariableNames', {'Description', 'Cost [USD]'});

% Display the table
fprintf('Table 6 (in %s): Total Cost Analysis (Fuel, Electricity, Opeation, Capital)\n',name_period_analyzed)
disp(Cost_Table);
end

