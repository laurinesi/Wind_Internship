% Block Maxima = ne prendre que la valeur maximale de chaque année


function [max_values] = BM_select(dataset)
    % Check if dataset is a table
    if ~istable(dataset)
        error('Input dataset must be a table.');
    end
    
    % Extract years and corresponding data from the table
    years = unique(dataset{:,1});  % Assuming the first column is 'Year'
    max_values = zeros(length(years), 2);  % Preallocate for speed
    
    for i = 1:length(years)
        year = years(i);
        
        % Extract data for the current year
        year_data = dataset{dataset{:,1} == year, 2};  % Assuming the second column contains the data
        
        % Check if year_data is empty to avoid errors
        if isempty(year_data)
            error('No data found for year: %d', year);
        end
        
        % Find the maximum value for the current year
        max_value = max(year_data);
        max_values(i, :) = [year, max_value];  % Store the year and max value
    end

    % Create a table with appropriate column names
    FF_BM = array2table(max_values, 'VariableNames', {'Year', 'FF'});

    % Save the table as a CSV file
    writetable(FF_BM, 'Schiphol_BM.csv');
    
    % Save the table as a text file
    writetable(FF_BM, 'Schiphol_BM.txt', 'Delimiter', '\t');

end


