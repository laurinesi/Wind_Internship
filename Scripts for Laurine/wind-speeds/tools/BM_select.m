% Block Maxima function
%
% max_values : table [years max_ws]

function [max_values] = BM_select(dataset)

% Check if dataset is a table
if ~istable(dataset)
    error('Input dataset must be a table.');
end

% Extract years and corresponding data from the table
years = unique(dataset.Year);
max_values = zeros(length(years), 2);

for i = 1:length(years)
    year = years(i);
    
    % Extract data for the current year
    idx = dataset.Year == year;
    year_data = dataset.F010(idx);
    
    % Find the maximum value for the current year
    max_value = max(year_data);
    max_values(i, :) = [year, max_value];  % Store the year and max value
end

% Create a table
FF_BM = array2table(max_values, 'VariableNames', {'Year', 'FF'});

% % Save the table as a CSV file
% writetable(FF_BM, 'Schiphol_BM.csv');

% Save the table as a text file
writetable(FF_BM, 'windspeeds_BM.txt', 'Delimiter', '\t');

end