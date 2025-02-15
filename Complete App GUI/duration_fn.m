function [days_per_month, hours_per_month, seconds_per_month, total_days, total_hours, total_seconds]=duration_fn(months_design)
    % Default Parameters 
    days_per_month_full = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % Full year
    months_full = {'January', 'February', 'March', 'April', ...
                   'May', 'June', 'July', 'August', ...
                   'September', 'October', 'November', 'December'};
    
    % Match the days_per_month for the selected winter months
    days_per_month = arrayfun(@(m) days_per_month_full(strcmp(months_full, m)), months_design);
    hours_per_month = days_per_month * 24; 
    seconds_per_month = hours_per_month * 3600; 
    % Calculate Total Days and Hours
    total_days = sum(days_per_month); % Total days in the winter season
    total_hours = total_days * 24; % Total hours in the winter season
    total_seconds = total_hours * 3600; % Total time in seconds
    
    fprintf('\nWinter Season Duration:\n');
    fprintf('Months: %s\n', strjoin(months_design, ', '));
    fprintf('Total Days: %d\n', total_days);
    fprintf('Total Hours: %d\n', total_hours);
    
end