function HospitalAirDistributionApp

    % Close any existing instances of the GUI
    existingApp = findall(0, 'Type', 'figure', 'Name', 'Hospital Air Distribution');
    if ~isempty(existingApp)
        delete(existingApp);
    end

    % Create the GUI figure
    app = uifigure('Name', 'Hospital Air Distribution', 'Position', [100, 100, 1000, 800]);

    % Dropdown for selecting design option
    uilabel(app, 'Position', [20, 750, 150, 20], 'Text', 'Select Design Option:');
    designOptionDropDown = uidropdown(app, ...
        'Position', [180, 750, 100, 22], ...
        'Items', {'Option 1 (Direct Heating by Furnace)', 'Option 2 (CHEHP)', 'Option 3 (CHEHP + HEngine Heat Recovery)', 'Option 4 (CHEHP + Heat Recovery from HEngine and Air Discharge))'}, ...
        'Value', 'Option 1 (Direct Heating by Furnace)', ...
        'ValueChangedFcn', @(~, ~)updateSchematic());

    % Input fields for user parameters
    uilabel(app, 'Position', [20, 700, 200, 20], 'Text', 'Ambient Temperature [C] (T0_list):');
    T0Field = uieditfield(app, 'text', 'Position', [230, 700, 200, 22], 'Value', '-10:4:10');

    uilabel(app, 'Position', [20, 660, 200, 20], 'Text', 'Ambient Pressure [bar] (P0):');
    P0Field = uieditfield(app, 'numeric', 'Position', [230, 660, 200, 22], 'Value', 1);

    uilabel(app, 'Position', [20, 620, 200, 20], 'Text', 'Design Temperature [C] (T3):');
    T3Field = uieditfield(app, 'numeric', 'Position', [230, 620, 200, 22], 'Value', 36);

    uilabel(app, 'Position', [20, 580, 200, 20], 'Text', 'Design Pressure [bar] (P3):');
    P3Field = uieditfield(app, 'numeric', 'Position', [230, 580, 200, 22], 'Value', 1.1);

    uilabel(app, 'Position', [20, 540, 200, 20], 'Text', 'Discharge Temperature [C] (T4):');
    T4Field = uieditfield(app, 'numeric', 'Position', [230, 540, 200, 22], 'Value', 24);

    uilabel(app, 'Position', [20, 500, 200, 20], 'Text', 'Discharge Pressure [bar] (P4):');
    P4Field = uieditfield(app, 'numeric', 'Position', [230, 500, 200, 22], 'Value', 1);

    % Button to run the selected design option
    runButton = uibutton(app, 'push', ...
        'Text', 'Run', ...
        'Position', [20, 460, 100, 30], ...
        'ButtonPushedFcn', @(~, ~)runDesignOption());

    % Table 1: Cost Table
    uilabel(app, 'Position', [500, 730, 250, 20], 'Text', 'Cost Summary', 'FontWeight', 'bold');
    resultsTable1 = uitable(app, ...
        'Position', [500, 500, 450, 200], ...
        'Data', {}, ...
        'ColumnName', {}, ...
        'RowName', {}, ...
        'Tag', 'Cost Table');

    % Table 2: Carbon Emissions Table
    uilabel(app, 'Position', [500, 470, 250, 20], 'Text', 'Carbon Emissions Summary', 'FontWeight', 'bold');
    resultsTable2 = uitable(app, ...
        'Position', [500, 350, 450, 100], ...
        'Data', {}, ...
        'ColumnName', {}, ...
        'RowName', {}, ...
        'Tag', 'Carbon Emissions Table');

    % UIAxes for the bar chart
    %[500, 260, 450, 20]
    ax1 = uiaxes(app, 'Position', [20, 20, 450, 300]); % Lowered chart
    title(ax1, 'Effect of Ambient Temperature Changes on Exergy Efficiency [%]');
    xlabel(ax1, 'Ambient Temperature (°C)');
    ylabel(ax1, 'Exergy Efficiency (%)');

    % Schematic Image Title
    uilabel(app, 'Position', [500, 260, 450, 20], ...
        'Text', 'Design Schematic Illustration', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % UIImage for the schematic diagram
    schematicImage = uiimage(app, ...
        'Position', [500, 50, 450, 180], ... % Enlarged image
        'ImageSource', 'default.png'); % Default schematic for Option 1

    % Callback function to update the schematic based on the selected option
    function updateSchematic()
        selectedOption = designOptionDropDown.Value;
        switch selectedOption
            case 'Option 1 (Direct Heating by Furnace)'
                schematicImage.ImageSource = 'Option1.png';
            case 'Option 2 (CHEHP)'
                schematicImage.ImageSource = 'Option2.png';
            case 'Option 3 (CHEHP + HEngine Heat Recovery)'
                schematicImage.ImageSource = 'Option3.png';
            case 'Option 4 (CHEHP + Heat Recovery from HEngine and Air Discharge))'
                schematicImage.ImageSource = 'Option4.png';
        end
    end

    % Callback function for the Run button
    function runDesignOption()
    % Get selected design option
    selectedOption = designOptionDropDown.Value;

    % Get user inputs
    try
        T0_list = eval(T0Field.Value); % Accepts a scalar or vector
        P0 = P0Field.Value;
        T3 = T3Field.Value;
        P3 = P3Field.Value;
        T4 = T4Field.Value;
        P4 = P4Field.Value;

        % Validate inputs
        if isempty(T0Field.Value) || isnan(P0) || isnan(T3) || ...
           isnan(P3) || isnan(T4) || isnan(P4)
            error('All fields must be filled with valid numeric values!');
        end

        % Execute the selected option and handle errors from relied functions
        try
            switch selectedOption
                case 'Option 1 (Direct Heating by Furnace)'
                    [Cost_Table, Carbon_Emissions_Table_winter, plotData] = Option1(T0_list, P0, T3, P3, T4, P4);
                case 'Option 2 (CHEHP)'
                    [Cost_Table, Carbon_Emissions_Table_winter, plotData] = Option2(T0_list, P0, T3, P3, T4, P4); % Placeholder
                case 'Option 3 (CHEHP + HEngine Heat Recovery)'
                    [Cost_Table, Carbon_Emissions_Table_winter, plotData] = Option3(T0_list, P0, T3, P3, T4, P4); % Placeholder
                case 'Option 4 (CHEHP + Heat Recovery from HEngine and Air Discharge))'
                    [Cost_Table, Carbon_Emissions_Table_winter, plotData] = Option4(T0_list, P0, T3, P3, T4, P4); % Placeholder
            end
        catch ME
            % Catch and display errors from relied functions
            uialert(app, ME.message, 'Simulation Error');
            return; % Stop further execution
        end

        % Display tables in the GUI
        resultsTable1.Data = table2cell(Cost_Table); % Populate table 1
        resultsTable1.ColumnName = Cost_Table.Properties.VariableNames; % Set headers for table 1
        resultsTable2.Data = table2cell(Carbon_Emissions_Table_winter); % Populate table 2
        resultsTable2.ColumnName = Carbon_Emissions_Table_winter.Properties.VariableNames; % Set headers for table 2

        % Plot the bar chart in the GUI axes
        cla(ax1); % Clear previous plot
        bar(ax1, plotData.x, plotData.y, 'grouped'); % Grouped bar chart
        legend(ax1, plotData.legend, 'Location', 'Best');
        title(ax1, 'Effect of Ambient Temperature Changes on Exergy Efficiency [%]');
        xlabel(ax1, 'Ambient Temperature (°C)');
        ylabel(ax1, 'Exergy Efficiency (%)');

    catch ME
        % Handle general errors (e.g., invalid inputs)
        uialert(app, ME.message, 'Input Error');
    end
end
end
